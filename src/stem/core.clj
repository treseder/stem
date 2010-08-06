(ns stem.core
  (:use [clojure.pprint] [clojure set])
  (:require [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as util]
            [stem.messages :as m]
            [stem.mle :as mle]
            [clojure.string :as str]
            [clojure.contrib.combinatorics :as comb])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml])
  (:gen-class))

(def *stem-version* 2.0)

;; if this value is true, then a caught exception will exit the
;; system.  If in dev, we don't want it to exit since it destroys the
;; REPL session.  This should be set to true before generating the
;; uber jar
(def in-production? false)


(defmacro with-exc [body message]
  `(try
    ~body
    (catch Exception e# (util/abort ~message e# ~in-production?))))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions to build up species and lineage maps and indexes into
;; maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn parse-settings-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files = a map - keys are filenames, values are the constant
  :species = a map - keys are specie names, values are comma separated lineage names
  :properties = a map with user configurable properties"
  ([] (parse-settings-file (util/get-settings-filename)))
  ([file-name]
     (let [f (if (nil? file-name) (util/get-settings-filename) file-name)
           yaml (Yaml.)
           yaml-map (with-exc (.load yaml (slurp f)) (m/e-strs :yaml))]
       (zipmap (map keyword (.keySet yaml-map))
               (map #(into {} %) (.values yaml-map))))))

(defn build-lin-to-spec-map
  "From the generic property map, builds the lineages to species map"
  [props]
  (reduce
   (fn [m [k v]]
     (let [lins  (str/split (util/remove-whitespace v) #",")]
       (merge m (zipmap lins (repeat k)))))
   {} (:species props)))

(defn build-spec-to-lin-map
  [props]
  (reduce
   (fn [m [k v]]
     (let [lins  (str/split (util/remove-whitespace v) #",")]
       (merge m {k (set lins)})))
   {} (:species props)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;functions for building the coalescent matrix ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn get-sorted-name-by-index [idx-map]
  (sort-by #(idx-map %) (keys idx-map)))

(defn create-name-to-index-map
  "The lineages names are strings that need to have indexes into the matrix.
  Returns a map of the n lineage names to indexes 0...n"
  [l-set]
  (zipmap l-set (range (count l-set))))

(defn name->indexes
  "Takes a seq of lineage names and turns them into the indexes of the matrix"
  [lins index-map]
  (map index-map lins))

(defn fill-matrix-for-lineage
  [i-idx j-idxs matrix c-time]
  (doseq [j-idx j-idxs]
    (if (> i-idx j-idx) ;; matrix is symmetric - only fill half
      (util/aset! matrix j-idx i-idx c-time)
      (util/aset! matrix i-idx j-idx c-time))))

(defn check-lin [idx name]
  (if (nil? idx) (util/abort (str (m/e-strs :missing-lin) name) nil in-production?))
  idx)

(defn rec-fill-time-matrix
  [node matrix lin-index parent-desc-set parent-c-time]
  (let [[{name :name  desc-set :desc c-time :c-time} l-node r-node] node
        ;; if leaf consider self the only desc
        checked-desc-set (if (empty? desc-set) #{name} desc-set)
        diff-set (difference parent-desc-set checked-desc-set)
        js (name->indexes diff-set lin-index)]
    (doseq [i checked-desc-set]
      (fill-matrix-for-lineage (check-lin (lin-index i) i) js matrix parent-c-time))
    (if-not (newick/is-leaf? node)
      (do (rec-fill-time-matrix l-node matrix lin-index desc-set c-time)
          (rec-fill-time-matrix r-node matrix lin-index desc-set c-time)))))

(defn fill-time-matrix-for-tree
  "Returns a matrix filled with coalescent times for each of the lineages.
  The matrix is upper triangular."
  [tree matrix lin-index]
  (rec-fill-time-matrix tree matrix lin-index (:desc (first tree)) 0)
  matrix)

(defn gene-trees-to-matrices
  "Takes a seq of gene-trees and returns a seq of matrices"
  [trees lin-indexes]
  (let [len (count lin-indexes)]
    (map #(fill-time-matrix-for-tree (:vec-tree %)
                                     (util/make-stem-array len len)
                                     lin-indexes)
         trees)))

(def not-zero? (complement zero?))

(defn reduce-matrices
  "Assumes matrices are pxp, upper-triangular.  Checks each cell in each
  matrix and stores the lesser in a.  Returns a."
  [a b]
  (let [p (count a)]
    (dorun
     (for [i (range p) j (range (+ i 1) p)]
       (let [a-val (util/aget! a i j)
             b-val (util/aget! b i j)]
         (if (or (zero? a-val) (and (< b-val a-val) (not-zero? b-val))  )
           (util/aset! a i j b-val)))))
    a))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions for generating species matrix;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn fill-matrix-for-species [i j c-time spec-to-index s-matrix]
  (let [current-val (if (> i j) (util/aget! s-matrix j i)  (util/aget! s-matrix i j) )]
    (if  (or (< c-time current-val) (= current-val 0))
      (if (< i j)
        (util/aset! s-matrix i j c-time)
        (util/aset! s-matrix j i c-time)))))

(defn to-spec-matrix
  [l-matrix {:keys [lin-to-index spec-set spec-to-index lin-to-spec]}]
  (let [s-count (count spec-set)
        s-matrix (util/make-stem-array s-count s-count)
        l-count (count lin-to-spec)
        index-to-lin (zipmap (vals lin-to-index) (keys lin-to-index))]
    ;; for loop only traverses the upper right triangle of the matrix
    (dorun (for [i (range l-count) j (range (+ i 1) l-count)]
       (let [i-lin (index-to-lin i)
             j-lin (index-to-lin j)
             i-spec (spec-to-index (lin-to-spec i-lin))
             j-spec (spec-to-index (lin-to-spec j-lin))]
         (if-not (= i-spec j-spec)
           (fill-matrix-for-species i-spec j-spec (util/aget! l-matrix i j)
                                    spec-to-index s-matrix)))))
    s-matrix))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions for generating newick str from species matrix;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn sort-vec-of-vec [coll]
  (sort #(compare (first %1) (first %2)) coll))

(defn get-list-permutations
  "The list of nodes used to build the ml species tree contains ties in how
  the tree should be built, with some nodes having the same coalescent time.
  Find all permutations of this list, thus finding all ml trees"
  [lst]
  (apply comb/cartesian-product (vals (group-by first lst))))

(defn matrix->sorted-list
  "Returns a list sorted by elements of the matrix.
  Assumes matrix is upper triangular where diag elements are zero"
  [m idx-to-name]
  (let [size (count m)
        v (transient [])]
    (dorun
     (for [i (range size) j (range (+ i 1) size)]
       (let [val (util/aget! m i j)]
         (conj! v [val (idx-to-name i) (idx-to-name j)]))))
    (sort #(compare (first %1) (first %2)) (persistent! v))))

(defn partial-find
  "Each key in m is a set of keys. Return value of map entry where k
  is a member of such set."
  [k m]
  (first (filter (fn [[ks v]] (contains? ks k)) m)))

(defn nil->node [n v]
  (if (vector? v) v [n]))

(defn nil->set [n aset]
  (if (nil? aset) #{n} aset))

(defn find-and-merge-nodes
  "Given a node and a map of nodes merge the node with its appropriate
  ancestor.
  node-map is of the form {#{s1 s2} [.55 [s1][s2]]}
  node is a vector of size 3 that contains the coalescent time of the two species:
  [0.55 s1 s2].  Ancestors are all internal nodes, with the species being leaves.
  of the tree."
  [node-map node]
  (let [[time l-name r-name] node
        [l-descs l-tree] (partial-find l-name node-map)
        [r-descs r-tree] (partial-find r-name node-map)]
    (if-not (or (contains? l-descs r-name)  (contains? r-descs l-name))
      (let [comb-tree [time (nil->node l-name l-tree) (nil->node r-name r-tree)]
            new-descs (union (nil->set l-name l-descs) (nil->set r-name r-descs))]
           (->
            node-map
            (dissoc l-descs r-descs)
            (assoc new-descs comb-tree)))
      node-map)))

(defn tree-from-seq [s]
  (reduce #(find-and-merge-nodes %1 %2) {} s))

(defmacro with-eval-and-out [form res-fn e-message & args]
  `(try
    (let [res# ~form]
      (~res-fn res# ~@args)
      res#)
    (catch Exception e#
      (util/abort ~e-message e# ~in-production?))))

(defprotocol JobProtocol
  (pre-run-check [job])
  (print-job [job])
  (run [job])
  (print-results [job])
  (print-results-to-file [job filename]))

(defrecord MLEJob [props env gene-trees results]
  JobProtocol
  (pre-run-check
   [job]
   (let [{:keys [props env gene-trees]} job]
     (doseq [k [:theta :lin-set :spec-set]]
       (util/abort-if-empty (env k) (m/e-strs k) in-production?)))
   job)

  (print-job
   [job]
   (m/print-job job)
   job)
  
  (run
   [job]
   (let [{:keys [props env gene-trees]} job
         gene-matrices  (gene-trees-to-matrices gene-trees (env :lin-to-index))
         min-gene-matrix (reduce reduce-matrices (util/make-stem-array (env :mat-size)) gene-matrices)
         spec-matrix (to-spec-matrix min-gene-matrix env)
         spec-lst (matrix->sorted-list spec-matrix (env :index-to-spec))
         lst-of-perm (get-list-permutations spec-lst)
         [tree-set tree] (first (tree-from-seq spec-lst))
         species-newick (newick/tree->newick-str tree)
         species-vec-tree (newick/build-tree-from-newick-str species-newick 1.0 1.0)
         mle (mle/calc-mle gene-trees species-vec-tree (env :spec-to-lin) (env :theta))
         res {:tied-trees lst-of-perm, :species-tree species-newick, :mle mle}]
     (assoc job :results res)))
  
  (print-results
   [job]
   (m/print-job-results (:results job))
   job)
  
  (print-results-to-file
   [job filename]
   (util/write-to-file filename ((job :results) :species-tree))
   job))

(defrecord SearchJob [props env gene-trees results])
(defrecord UserTreeJob [props env gene-trees results])

(defn create-env
  "This function returns a hashmap of all of the various data-structures that are necessary
  to generate species trees and their likelihoods"
  [settings-map]
  (let [lin-to-spec (build-lin-to-spec-map settings-map)
        spec-to-lin (build-spec-to-lin-map settings-map)
        lin-set (set (keys lin-to-spec))
        spec-set (set (vals lin-to-spec))
        spec-to-index (create-name-to-index-map spec-set)
        index-to-spec (zipmap (vals spec-to-index) (keys spec-to-index))        
        lin-to-index (create-name-to-index-map lin-set)]
    {:theta ((:properties settings-map) "theta")
     :lin-to-spec lin-to-spec
     :spec-to-lin spec-to-lin
     :lin-set lin-set
     :spec-set spec-set
     :spec-to-index spec-to-index
     :index-to-spec index-to-spec
     :lin-to-index lin-to-index
     :mat-size (count lin-set)}))

(defn create-job [& file]
  (let [{:keys [properties files species] :as s}
        (parse-settings-file (first file))
        env (create-env s)
        gene-trees (g-tree/get-gene-trees files (env :theta))]
    (case (properties "run")
          "0" (UserTreeJob. properties env gene-trees nil)
          "2" (SearchJob. properties env gene-trees nil)
          (MLEJob. properties env gene-trees nil))))

(defn -main
  "Entry point for STEM 2.0. Throughout the code, lin refers to lineages, and spec refers
  to species."
  [& args]
  (m/header-message *stem-version*)
  (-> (create-job) (pre-run-check) (print-job) (run) (print-results))
  (if in-production? (System/exit 0)))
