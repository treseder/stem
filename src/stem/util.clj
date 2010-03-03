(ns stem.util
  (:require [clojure.contrib.str-utils2 :as s])
  (:import [java.io File StringReader BufferedReader]))

(defn abort [message e & exit]
  (println message)
  (if-not (nil? e)
    (println (.getMessage e)))
  (if exit (System/exit 1)))


(defmulti to-float class)

(defmethod to-float clojure.lang.Keyword [k-word]
  (Double/parseDouble (name k-word)))

(defmethod to-float java.lang.String [s]
  (Double/parseDouble s))


(defn get-file
  [f-name]
  (try
   (File. f-name)
   (catch Exception _
     (abort (str "Couldn't find the file: " f-name)))))


(defn remove-whitespace [s]
  (s/replace s #"\s" ""))
