(ns stem.util
  (:require [clojure.contrib.str-utils2 :as s])
  (:import [java.io File]))

(defn abort [message & exit?]
  (println message)
  (if exit? (System/exit 1)))


(defn get-file
  [f-name]
  (try
   (File. f-name)
   (catch Exception _
     (abort (str "Couldn't find the file: " f-name)))))


(defn remove-whitespace [s]
  (s/replace s #"\s" ""))

