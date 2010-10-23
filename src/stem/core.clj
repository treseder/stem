(ns stem.core
  (:use [stem.constants])
  (:require [stem.util :as util]
            [stem.messages :as m]
            [stem.job :as j])
  (:gen-class))

(defn -main
  "Entry point for STEM 2.0. Throughout the code, lin refers to lineages, and spec refers
  to species."
  [& args]
  (m/header-message *stem-version*)
  (-> (j/create-job) (j/pre-run-check) (j/print-job) (j/run) (j/print-results)
      (j/print-results-to-file))
  (m/print-done)
  (if *in-production* (System/exit 0)))
