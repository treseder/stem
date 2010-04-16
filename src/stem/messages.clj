(ns stem.messages
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.util :as util]))

(def e-strs {:yaml "An error occurred parsing the settings file."})

(defn yaml-message [prop-map]
  (println "Successfully parsed settings.yaml config file")
  (flush))

(defn theta-message [theta]
  (println "Using theta = " theta) (println) (flush))

(defn lin-set-message [s]
  (println "There are " (count s) " lineages:")
  (println (apply str (interpose "," s)))
  (println) (flush))

(defn spec-set-message [s]
  (println "There are" (count s) "species:")
  (println (apply str (interpose "," s)))
  (println) (flush))
