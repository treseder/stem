(ns stem.messages
  (:use [stem.constants])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.util :as util]))

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

(defn print-lik-job-results [{:keys [tied-trees species-tree mle]}]
  (let [tt (count tied-trees)]
    (println "\n\n****************Results*****************")
    (println (str "\nLikelihood Species Tree (newick format):\n\n" species-tree))
    (println-if (> tt 1) (str "\nNote: There were " tt " trees that have the same likelihood estimate as the tree above.  These trees will be output to the 'stem-tree.tre' file."))
    (println (str "\nLikelihood estimate for tree: " mle))))

(defn header-message [version]
  (println      "***************************************")
  (println (str "**        Welcome to STEM " version "        **"))
  (println      "***************************************\n")
  (flush))

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
  (println (str "Using burnin: " (get props "burnin" *burnin-default*)))
  (println (str "Using bound_total_iters: " (get props "bound_total_iter" *bound-total-iter-default*)))
  (println (str "Using num_saved_trees:: " (get props "num_saved_trees" *num-saved-trees-default*)))
  (println "\nBeginning search now (this could take a while)...\n"))

(defn print-user-job [{:keys [props env gene-trees]}]
  (println "Species tree read from 'user.tre':")
  (println)
  (println (env :species-newick)))

(defn print-search-results [{:keys [best-trees]}]
  (println "Search completed.")
  (println "Here are the results (also written to file 'search.trees':")
  (doseq [[lik n-str] best-trees] (println (str "[" lik "] " n-str))))

(defn print-user-results [{:keys [likelihood]}]
  (println (str "\nThe likelihood for the tree is: " likelihood)))

(defn yaml-message [prop-map]
  (println "Successfully parsed the settings file")
  (flush))

(defn theta-message [theta]
  (println "Using theta =" theta) (println) (flush))

(defn lin-set-message [s]
  (println "There are" (count s) "lineages:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

(defn spec-set-message [s]
  (println "There are" (count s) "species:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

(defn spec-newick-message [s num-trees]
  (println "The Maximum Likelihood tree (written to file 'mle.tre') is:")
  (println)
  (println s)
  (println-if (> num-trees 1) (str "\nThere were also " num-trees " trees that have the same likelihood estimate as the tree above.\nThese trees are also output to 'mle.tre'."))
  (flush))

(defn mle-message [mle]
  (println (str "The likelihood estimate is: " mle))
  (flush))
