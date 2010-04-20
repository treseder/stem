(ns stem.core
  (:use [clojure.contrib pprint] [clojure set])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as util]
            [stem.messages :as m])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml])
  (:gen-class))

;; if this value is true, then a caught exception will exit the
;; system.  If in dev, we don't want it to exit since it destroys the
;; REPL session.  This should be set to true before generating the
;; uber jar
(def in-production? false)

(defn parse-yaml-config-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files = a map - keys are filenames, values are the constant
  :species = a map - keys are specie names, values are comma separated lineage names
  :properties = a map with user configurable properties"
  [file-name]
  (let [yaml (Yaml.)
        yaml-map (.load yaml (slurp file-name))]
    (zipmap (map keyword (.keySet yaml-map))
            (map #(into {} %) (.values yaml-map)))))

(defn get-lineages-to-spec-map
  "From the generic property map, builds the lineages to species map"
  [map]
  (let [s-map (:species map)]
    (reduce
     (fn [m [k v]]
       (let [lineages  (s/split (util/remove-whitespace v) #",")]
         (merge m (zipmap lineages (repeat k)))))
     {} s-map)))

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
  [i js matrix c-time]
  (doseq [j js]
    (if (> i j) ;; matrix is symmetric - only fill half
      (util/aset! matrix j i c-time)
      (util/aset! matrix i j c-time))))

(defn rec-fill-time-matrix
  [node matrix lin-index parent-desc-set parent-c-time]
  (let [[{name :name  desc-set :desc c-time :c-time} l-node r-node] node
        ;; if leaf consider self the only desc
        checked-set (if (empty? desc-set) #{name} desc-set)
        diff-set (difference parent-desc-set checked-set)
        js (name->indexes diff-set lin-index)]
    (doseq [i checked-set] (fill-matrix-for-lineage (lin-index i) js matrix parent-c-time))
    (if-not (newick/is-leaf? node)
      (do (rec-fill-time-matrix l-node matrix lin-index desc-set c-time)
          (rec-fill-time-matrix r-node matrix lin-index desc-set c-time)))))

(defn fill-time-matrix-for-tree
  "Returns a matrix filled with coalescent times for each of the lineages.
  The matrix is upper triangular."
  [tree matrix lin-index]
  (rec-fill-time-matrix tree matrix lin-index ((first tree) :desc) 0)
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

(defn to-species-matrix [l-matrix lin-to-index spec-set spec-to-index lin-to-spec]
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

(defmacro with-exc-and-out [form res-f e-message & args]
  `(try
    (let [res# ~form]
      (if (nil? ~args)
        (~res-f res#)
        (~res-f res# ~@args))
      res#)
    (catch Exception e#
      (util/abort ~e-message e# ~in-production?))))

(defn -main
  "Entry point for STEM 2.0. Throughout the code, lin refers to lineages, and spec refers
  to species."
  [& args]
  (let [prop-map (with-exc-and-out (parse-yaml-config-file "settings.yaml") m/yaml-message (m/e-strs :yaml))
        theta (with-exc-and-out ((:properties prop-map) "theta") m/theta-message "")
        lin-to-spec (get-lineages-to-spec-map prop-map)
        lin-set (with-exc-and-out (set (keys lin-to-spec)) m/lin-set-message "")
        spec-set (with-exc-and-out (set (vals lin-to-spec)) m/spec-set-message "")
        spec-to-index (create-name-to-index-map spec-set)
        index-to-spec (zipmap (vals spec-to-index) (keys spec-to-index))        
        lin-to-index (create-name-to-index-map lin-set)
        gene-trees (g-tree/get-gene-trees (:files prop-map) theta)
        matrices  (gene-trees-to-matrices gene-trees lin-to-index)
        m-size (count lin-set)
        least-matrix (reduce reduce-matrices (util/make-stem-array m-size m-size) matrices)
        species-matrix (to-species-matrix least-matrix lin-to-index spec-set spec-to-index lin-to-spec)
        lst (matrix->sorted-list species-matrix index-to-spec)
        [tree-set tree] (first (tree-from-seq lst))
        species-newick (with-exc-and-out (newick/tree->newick-str tree) m/spec-newick-message "")]
    (println lst)
    (util/write-to-file "mletree.tre" species-newick))
  (if in-production? (System/exit 0)))
