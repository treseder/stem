(ns stem.util
  (:require [clojure.contrib.str-utils2 :as s])
  (:import [java.io File StringReader BufferedReader]))

(defn abort [message e & exit]
  (println message)
  (if-not (nil? e)
    (println (.getMessage e)))
  (if exit (System/exit 1)))


(defmulti to-double class)

(defmethod to-double clojure.lang.Keyword [k-word]
  (Double/parseDouble (name k-word)))

(defmethod to-double java.lang.String [s]
  (Double/parseDouble s))

(defmethod to-double java.lang.Integer [i]
  (double i))

(defn get-file
  [f-name]
  (try
   (File. f-name)
   (catch Exception _
     (abort (str "Couldn't find the file: " f-name)))))

(defn remove-whitespace [s]
  (s/replace s #"\s" ""))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; multi-dim array functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn array? [x] (-> x class .isArray))

(defn see-array [x] (if (array? x) (println (map see-array x)) x))

(defn array->str [a]
 (str "[" (reduce #(str %1 " " %2) a) "]"))

(defn print-array [a]
  "Prints a 1D or 2D array"
  (println (if (> (count a) 1)
     (reduce #(str %1 "\n" (array->str %2)) "" a)
     (array->str a))))

(defn make-stem-array [rows cols]
  (make-array Double/TYPE rows cols))

(defmacro aget!
  ([array y]      `(aget ~(vary-meta array assoc :tag 'doubles) ~y))
  ([array x & ys] `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~@ys)]
                     (aget! a# ~x))))

(defmacro aset! [array x y v]
  (let [nested-array `(aget ~(vary-meta array assoc :tag 'objects) ~y)
        a-sym         (with-meta (gensym "a") {:tag 'doubles})]
    `(let [~a-sym ~nested-array]
       (aset ~a-sym ~x (double ~v)))))
