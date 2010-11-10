(ns stem.job
  (:use [stem.constants])
  (:require [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as u]
            [stem.messages :as m]
            [stem.lik :as lik]
            [clojure.string :as str]
            [stem.lik-tree :as lt]
            [stem.search :as s])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml]))


(defprotocol JobProtocol
  (pre-run-check [job])
  (print-job [job])
  (run [job])
  (search [job])
  (print-results [job])
  (print-results-to-file [job]))


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
         mle-newick (newick/tree->newick-str (first tied-trees))
         mle-vec-tree (newick/build-tree-from-newick-str mle-newick 1.0 1.0)
         mle (lik/calc-mle gene-trees mle-vec-tree (env :spec-to-lin) (env :theta))
         user-vec-trees (map #(newick/build-tree-from-newick-str % 1.0 1.0) (env :user-trees))
         user-liks (lik/calc-mles gene-trees user-vec-trees (env :spec-to-lin) (env :theta))
         res {:user-liks user-liks :mle mle :mle-newick mle-newick}]
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
         mle (lik/calc-mle gene-trees species-vec-tree (env :spec-to-lin) (env :theta))
         res {:tied-trees tied-trees, :species-matrix spec-matrix, :species-tree species-newick, :mle mle}]
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
          (map (fn [[lik vec-tree]] [lik (newick/vector-tree->newick-str vec-tree)]))
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
  :files = a map - keys are filenames, values are the constant
  :species = a map - keys are species names, values are comma separated lineage names
  :properties = a map with user configurable properties"
  ([] (parse-settings-file (u/get-settings-filename)))
  ([file-name]
     (let [f (if (nil? file-name) (u/get-settings-filename) file-name)
           yaml (Yaml.)
           yaml-map (u/with-exc (.load yaml (slurp f)) (m/e-strs :yaml))
           s-map (zipmap (map keyword (.keySet yaml-map))
                         (map #(into {} %) (.values yaml-map)))]
       (doseq [k [:properties :species]]
         (u/abort-if-empty (k s-map) (k m/e-strs)))
       s-map)))

(defn parse-user-tree-file [f]
  (-> (u/with-exc (u/read-file f) "You must specify a species tree in 'user.tre' to compute the likelihood.")
      (u/remove-whitespace)
      (str/split #";")))

(defn create-job [& file]
  (let [{:keys [properties files species] :as s} (parse-settings-file (first file))
        env (create-env s)
        gene-trees (g-tree/get-gene-trees files (env :theta))]
    (case (properties "run")
          0 (UserTreeJob. properties (assoc env :user-trees (parse-user-tree-file "user.tre")) gene-trees nil)
          2 (SearchJob. properties env gene-trees nil)
          (LikJob. properties env gene-trees nil))))
