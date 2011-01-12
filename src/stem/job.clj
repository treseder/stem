(ns stem.job
  (:use [stem.constants])
  (:require [stem.newick :as newick]
            [stem.util :as u]
            [stem.messages :as m]
            [stem.lik :as lik]
            [clojure.string :as str]
            [stem.lik-tree :as lt]
            [stem.search :as s]
            [stem.gene-tree :as g-tree]
            [stem.hybrid :as h])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml]))


(defprotocol JobProtocol
  (pre-run-check [job])
  (print-job [job])
  (run [job])
  (search [job])
  (print-results [job])
  (print-results-to-file [job]))

(defrecord HybridJob [props env gene-trees results]
  JobProtocol
  (pre-run-check
   [job]
   (doseq [k [:theta :lin-set :spec-set]]
     (u/abort-if-empty (env k) (m/e-strs k)))
   job)
  
  (print-job
   [job]
   (m/print-hyb-job job)
   job)
  
  (run
   [job]
   (let [{:keys [spec-matrix]} (lt/get-lik-tree-parts gene-trees env false)
         h-specs (str/split (u/remove-whitespace (get props "hybrid_species" )) #",")
         optim-c-time-fn (partial newick/optimized-c-time (env :spec-to-index) spec-matrix)
         optim-hybrid-tree (h/fix-tree-times
                            (newick/build-tree-from-newick-str
                             (env :hybrid-newick) 1.0 1.0 optim-c-time-fn)
                            spec-matrix
                            (env :spec-to-index)
                            (set h-specs))
         hybrid-trees (map #(h/fix-tree-times % spec-matrix (env :spec-to-index) (set h-specs))
                           (h/make-hybrid-trees-for-specs optim-hybrid-tree h-specs))
         hybrid-newicks (map newick/vector-tree->newick-str hybrid-trees)
         gammas (h/find-gammas gene-trees hybrid-trees (env :spec-to-lin) (env :theta))
         k (h/compute-k (count h-specs) (count (env :spec-to-index)))
         aic (h/compute-aic (first gammas) k)
         res {:hybrid-trees hybrid-trees, :species-matrix spec-matrix
              :gammas gammas :k k :aic aic :hybrid-newicks hybrid-newicks}]
     (assoc job :results res)))
  
  (print-results
   [job]
   (m/print-hyb-results results)
   job)
  
  (print-results-to-file
   [job]
   job))


(defrecord UserTreeJob [props env gene-trees results]
  JobProtocol

  (pre-run-check
   [job]
   (doseq [k [:theta :lin-set :spec-set]]
     (u/abort-if-empty (env k) (m/e-strs k)))
   job)
  
  (print-job
   [job]
   (m/print-lik-job job)
   (m/print-user-job job)
   job)
  
  (run
   [job]
   (let [{:keys [tied-trees spec-matrix]} (lt/get-lik-tree-parts gene-trees env)
         optim-c-time-fn (partial newick/optimized-c-time (env :spec-to-index) spec-matrix)
         optim-vec-trees (map #(s/fix-tree-times
                                (newick/build-tree-from-newick-str % 1.0 1.0 optim-c-time-fn)
                                spec-matrix
                                (env :spec-to-index))
                              (env :user-trees))
         optim-liks (lik/calc-liks gene-trees optim-vec-trees (env :spec-to-lin) (env :theta))
         optim-newicks (map newick/vector-tree->newick-str optim-vec-trees)
         user-vec-trees (map #(newick/build-tree-from-newick-str % 1.0 1.0) (env :user-trees))
         user-liks (lik/calc-liks gene-trees user-vec-trees (env :spec-to-lin) (env :theta))
         res {:user-liks user-liks :optim-trees optim-newicks :optim-liks optim-liks}]
     (assoc job :results res)))
  
  (print-results
   [job]
   (m/print-user-results (env :user-trees) results)
   job)
  
  (print-results-to-file
   [job]
   job))

(defrecord LikJob [props env gene-trees results]
  JobProtocol
  (pre-run-check
   [job]
   (doseq [k [:theta :lin-set :spec-set]]
     (u/abort-if-empty (env k) (m/e-strs k)))
   job)

  (print-job
   [job]
   (m/print-lik-job job)
   job)
  
  (run
   [job]
   (let [{:keys [min-gene-matrix tied-trees spec-matrix]} (lt/get-lik-tree-parts gene-trees env)
         species-newick (newick/tree->newick-str (first tied-trees))
         species-vec-tree (newick/build-tree-from-newick-str species-newick 1.0 1.0)
         mle (lik/calc-lik gene-trees species-vec-tree (env :spec-to-lin) (env :theta))
         res {:tied-trees tied-trees, :species-matrix spec-matrix,
              :species-tree species-newick, :mle mle}]
     (assoc job :results res)))
  
  (print-results
   [job]
   (m/print-lik-job-results results)
   job)
  
  (print-results-to-file
   [job]
   (let [f (get props "mle-filename" *mle-filename-default*)]
     (u/write-lines f (map
                       #(newick/tree->newick-str %)
                       (:tied-trees results)))
     (u/write-lines *bootstrap-filename* (map
                       #(newick/tree->newick-str % (str "#" (env :theta)))
                       (:tied-trees results))))))

(defrecord SearchJob [props env gene-trees results]
  JobProtocol
  (pre-run-check
   [job]
   job)
  
  (print-job
   [job]
   (m/print-search-job job)
   job)

  (run
   [job]
   (let [{:keys [tied-trees spec-matrix]} (lt/get-lik-tree-parts gene-trees env)
         species-newick (newick/tree->newick-str (first tied-trees))
         species-vec-tree (newick/build-tree-from-newick-str species-newick 1.0 1.0)]
     (->> (s/search-for-trees species-vec-tree gene-trees spec-matrix props env)
          (map (fn [sr] [(.lik sr) (newick/vector-tree->newick-str (.tree sr))]))
          (hash-map :best-trees)
          (assoc job :results))))
  
  (print-results
   [job]
   (m/print-search-results results)
   job)
  
  (print-results-to-file
   [job]
   (let [f (get props "search-filename" *search-filename-default*)]
     (u/write-lines f (map
                       (fn [[lik newick]] (str "[" (u/format-time lik) "]" newick))
                       (:best-trees results))))
   job))

(defn build-lin-to-spec-map
  "From the generic property map, builds the lineages to species map"
  [props]
  (reduce
   (fn [m [k v]]
     (let [lins  (str/split (u/remove-whitespace v) #",")]
       (merge m (zipmap lins (repeat (str k))))))
   {} (:species props)))

(defn build-spec-to-lin-map
  [props]
  (reduce
   (fn [m [k v]]
     (let [lins  (str/split (u/remove-whitespace v) #",")]
       (merge m {(str k) (set lins)})))
   {} (:species props)))

(defn create-name-to-index-map
  "The lineages names are strings that need to have indexes into the matrix.
  Returns a map of the n lineage names to indexes 0...n"
  [l-set]
  (zipmap l-set (range (count l-set))))

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

(defn parse-settings-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files = a map - keys are filenames, values are the rate constant
  :species = a map - keys are species names, values are comma separated lineage names
  :properties = a map with user configurable properties"
  ([] (parse-settings-file (u/get-settings-filename)))
  ([file-name]
     (let [f (if (nil? file-name) (u/get-settings-filename) file-name)
           yaml (Yaml.)
           yaml-map (u/with-exc
                      (.load yaml (str/replace (slurp f) "\t" "  ")) (m/e-strs :yaml))
           s-map (zipmap (map keyword (.keySet yaml-map))
                         (map #(into {} %) (.values yaml-map)))]
       (doseq [k [:properties :species]]
         (u/abort-if-empty (k s-map) (k m/e-strs)))
       s-map)))

(defn parse-user-tree-file [f]
  (-> (u/with-exc (u/read-file f) (str "You must specify a species tree in '" f "' to compute the likelihood."))
      (u/remove-whitespace)
      (str/split #";")))

(defn create-job [& file]
  (let [{:keys [properties files species] :as s} (parse-settings-file (first file))
        env (create-env s)
        gene-trees (g-tree/get-gene-trees files (env :theta))]
    (case (properties "run")
          0 (UserTreeJob. properties
                          (assoc env :user-trees (parse-user-tree-file "user.tre")) gene-trees nil)
          2 (SearchJob. properties env gene-trees nil)
          3 (HybridJob. properties
                        (assoc env :hybrid-newick (first (parse-user-tree-file "hybrid.tre")))
                        gene-trees nil)
          (LikJob. properties env gene-trees nil))))
