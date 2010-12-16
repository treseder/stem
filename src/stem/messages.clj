(ns stem.messages
  (:use [stem.constants])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.util :as util]
            [stem.newick :as n]))

(def e-strs {:yaml "An error occurred parsing the settings file."
             :theta "A value for theta must be present in the settings file."
             :lin-set "At least one species to lineage mapping must be present in the settings file"
             :spec-set "At least one species to lineage mapping must be present in the settings file"
             :missing-lin "A lineage encountered in the genetree file was not present in the settings file: "
             :properties "The properties section is missing from the settings file."
             :species "The species section is missing from the settings file."})

(defmacro println-if [test form]
  `(if ~test
    (println ~form)
    (flush)))
(defn header-message [version]
  (println      "***************************************")
  (println (str "**        Welcome to STEM " version "        **"))
  (println      "***************************************\n")
  (flush))

(defn print-done []
  (println "\n****************** Done ****************\n"))

(defn print-lik-job-results [{:keys [tied-trees species-tree species-matrix mle]}]
  (let [tt (count tied-trees)]
    (println "\n\n****************Results*****************\n")
    (println "D_AB Matrix:")
    (util/print-array species-matrix)
    (println (str "\nMaximum Likelihood Species Tree (Newick format):\n\n" species-tree))
    (println-if (> tt 1) (str "\nNOTE: There were at least" tt " trees that have the same log likelihood.  These trees will be output to the 'mle.tre' file."))
    (println (str "\nlog likelihood for tree: " mle))))

(defn print-hyb-job [{:keys [props env gene-trees]}]
  (println "The settings file was successfully parsed...\n")
  (println (str "Using theta = " (env :theta) "\n"))
  (println (str "The settings file contained " (count (env :spec-set)) " species and " (count (env :lin-set)) " lineages.\n"))
  (println "The species-to-lineage mappings are:\n")
  (doseq [[k v] (env :spec-to-lin)] (println (str k ": " (apply str (interpose ", " v)))))
  (println "\nHybridization params:\n")
  (println (str "Hybrid species: " (get props "hybrid_species")))
  (println (str "\nHybrid input tree:\n"))
  (println (env :hybrid-newick)))

(defn print-hyb-results
  [{:keys [species-matrix hybrid-trees gammas k aic]}]
  (println "\n\n****************Results*****************\n")
  (println "Hybrid resolution trees:\n")
  (doseq [t hybrid-trees] (println (n/vector-tree->newick-str t)))
  (println "\nD_AB Matrix:")
  (util/print-array species-matrix)  
  (println (str "\nMaximum log likelihood: " (first gammas) "\n"))
  (let [idv (map vector (iterate inc 1) (second gammas))]
    (doseq [[idx val] idv] (println (str "gamma" idx ": " val ))))
  (println (str "\nk: " k))
  (println (str "\nAIC: " aic)))

(defn print-lik-job [{:keys [props env gene-trees]}]
  (println "The settings file was successfully parsed...\n")
  (println (str "Using theta = " (env :theta) "\n"))
  (println (str "The settings file contained " (count (env :spec-set)) " species and " (count (env :lin-set)) " lineages.\n"))
  (println "The species-to-lineage mappings are:\n")
  (doseq [[k v] (env :spec-to-lin)] (println (str k ": " (apply str (interpose ", " v))))))

(defn print-search-job [{:keys [props env gene-trees]}]
  (println "The settings file was successfully parsed...\n\n")
  (println "Setting used for STEM search:")
  (println (str "Using beta: " (get props "beta" *beta-default*)))
  (println-if (get props "seed") (str "Using seed: " (get props "seed")))
  (println (str "Using burnin: " (get props "burnin" *burnin-default*)))
  (println (str "Using bound_total_iter: " (get props "bound_total_iter" *bound-total-iter-default*)))
  (println (str "Using num_saved_trees: " (get props "num_saved_trees" *num-saved-trees-default*)))
  (println "\nBeginning search now (this could take a while)...\n"))

(defn print-search-results [{:keys [best-trees]}]
  (println "Search completed.\n")
  (println "Here are the results (also written to file 'search.tre'):\n")
  (doseq [[lik n-str] best-trees] (println (str "[" (util/format-time lik) "] " n-str))))

(defn print-user-job [{:keys [props env gene-trees]}]
  (println (str "\nRead " (count (env :user-trees)) " species tree[s] from 'user.tre'")))

(defn print-user-results [user-newicks {:keys [user-liks optim-trees optim-liks]}]
  (println "\n\n****************Results*****************\n")
  (doseq [[tree lik] (partition 2 (interleave user-newicks user-liks))]
    (println "User tree: ")
    (println tree)
    (println (str "log likelihood for tree: " lik "\n")))
  (println "\n**************Optimized Trees************\n")
  (doseq [[tree lik] (partition 2 (interleave optim-trees optim-liks))]
    (println "Optimized user tree: ")
    (println tree)
    (println (str "log likelihood: " lik "\n"))))

(defn lin-set-message [s]
  (println "There are" (count s) "lineages:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

(defn spec-set-message [s]
  (println "There are" (count s) "species:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

