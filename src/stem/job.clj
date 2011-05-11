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
  (:import [java.io BufferedReader FileReader]))


(defprotocol JobProtocol
  (pre-run-check [job])
  (print-job [job])
  (run [job])
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
         spec->idx (env :spec-to-index)
         h-specs (str/split (u/remove-whitespace (get props "hybrid_species" )) #",")
         h-spec->idx (h/h-specs->gamma-ids h-specs)
         optim-c-time-fn (partial newick/optimized-c-time spec->idx  spec-matrix)
         user-parental-tree (newick/newick->tree (env :hybrid-newick) 1.0 1.0 optim-c-time-fn)
         parental-trees (map #(s/fix-tree-times % spec-matrix spec->idx)
                             (h/make-parental-trees user-parental-tree (set h-specs)))
         parental-data (h/parental-trees->parental-data parental-trees gene-trees
                                                        (env :spec-to-lin) (env :theta))
         
         ;; fixes parental tree lengths for hybridization
         gamma-top (h/gamma-topologys (set h-specs))
         hybrid-data (map #(h/gamma-topology->hybrid-tree-data
                            % parental-trees gene-trees spec-matrix h-spec->idx spec->idx
                            (env :spec-to-lin) (env :theta))
                          gamma-top)
         res {:species-matrix spec-matrix :hybrid-data hybrid-data :parental-data parental-data :h-spec->idx h-spec->idx}]
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
                                (newick/newick->tree % 1.0 1.0 optim-c-time-fn)
                                spec-matrix
                                (env :spec-to-index))
                              (env :user-trees))
         optim-liks (lik/calc-liks gene-trees optim-vec-trees (env :spec-to-lin) (env :theta))
         optim-newicks (map newick/vector-tree->newick-str optim-vec-trees)
         user-vec-trees (map #(newick/newick->tree % 1.0 1.0) (env :user-trees))
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
         species-vec-tree (newick/newick->tree species-newick 1.0 1.0)
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
         species-vec-tree (newick/newick->tree species-newick 1.0 1.0)]
     (->> (s/search-for-trees species-vec-tree gene-trees spec-matrix props env)
          (map (fn [[lik tree]] [lik (newick/vector-tree->newick-str tree)]))
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

(defn create-job [& file]
  (let [s (-> (u/parse-settings-file (first file) (m/e-strs :yaml))
              (u/check-settings-map m/e-strs))
        {:keys [properties files species]} s 
        env (create-env s)
        gene-trees (g-tree/get-gene-trees files (env :theta))]
    (case (properties "run")
          0 (UserTreeJob. properties
                          (assoc env :user-trees (u/parse-tree-file "user.tre")) gene-trees nil)
          2 (SearchJob. properties env gene-trees nil)
          3 (HybridJob. properties
                        (assoc env :hybrid-newick (first (u/parse-tree-file "hybrid.tre")))
                        gene-trees nil)
          (LikJob. properties env gene-trees nil))))

(defn run-job
  [file]
  (do (-> (create-job file)
          (pre-run-check)
          (print-job)
          (run)
          (print-results))
      (println "Finished")))
