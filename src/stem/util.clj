(ns stem.util
  (:use [clojure.pprint] [stem.constants])
  (:require [clojure.string :as str]
            [clojure.java.io :as io])
  (:import [java.io File StringReader BufferedReader]
           [java.util Random]))

(defmacro with-exc [body message]
  `(try
    ~body
    (catch Exception e# (stem.util/abort ~message e#))))

(defn rand-generator
  [& [seed]]
  (let [rg (if seed (Random. seed) (Random.))]
    (fn [& [type n]] (if (= type :int) (.nextInt rg n) (.nextDouble rg)))))


(def not-zero? (complement zero?))

(defn comma-sep [seq]
  (apply str (interpose "," seq)))

(defn println-flush [s]
  (println s)
  (flush))

(defn make-counter [init-val] 
  (let [c (atom init-val)] 
    {:next #(swap! c inc)
     :reset #(reset! c init-val)}))

(defn abort [message e]
  (println message)
  (if-not (nil? e)
    (println (.getMessage e)))
  (if *in-production* (System/exit 1)))

(defmulti abort-if-empty
  (fn [v m]
    (class v)))

(defmethod abort-if-empty java.lang.String
  [v message]
  (when (= v "") (abort message nil)))

(defmethod abort-if-empty java.lang.Number
  [v message])

(defmethod abort-if-empty :default
  [v message]
  (if (or (nil? v) (empty? v))
    (abort message nil)))

(defmulti to-double class)

(defmethod to-double clojure.lang.Keyword [k-word]
  (Double/parseDouble (name k-word)))

(defmethod to-double java.lang.String [s]
  (Double/parseDouble s))

(defmethod to-double java.lang.Integer [i]
  (double i))

(defmethod to-double java.lang.Double [i]
  i)

(defn read-file [f]
  (slurp f))

(defn write-to-file [f str]
  (spit f str))

(defn write-lines
  "Writes lines (a seq) to f, separated by newlines.  f is opened with
  writer, and automatically closed at the end of the sequence."
  [f lines]
  (with-open [writer (io/writer f)]
    (loop [lines lines]
      (when-let [line (first lines)]
        (.write writer (str line))
        (.newLine writer)
        (recur (rest lines))))))

(defn format-time [time]
  (cl-format nil "~,5f" time))

(defn get-file
  [f-name]
  (try
   (File. f-name)
   (catch Exception _
     (abort (str "Couldn't find the file: " f-name)))))

(defn contains-settings?
  [filename]
  (or (. filename equals "settings")
      (. filename equals "settings.txt")
      (. filename  equals "settings.yaml")))

(defn get-settings-filename
  "Searches the current directory for a file that is named
  'settings', or 'settings.txt', or 'settings.yaml', in that order"
  []
  (first (sort (filter contains-settings? (. (File. ".") list)))))


(defn remove-whitespace [s]
  (str/replace s #"\s" ""))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; multi-dim array functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn array? [x] (-> x class .isArray))

(defn see-array [x] (if (array? x) (map see-array x) x))

(defn array->str [a]
 (str "[" (reduce #(str %1 " " (format-time %2)) "" a) "]"))

(defn print-array [a]
  "Prints a 1D or 2D array"
  (println (if (> (count a) 1)
     (reduce #(str %1 "\n" (array->str %2)) "" a)
     (array->str a))))

(defn print-arrays [arrs]
  (doseq [a arrs] (print-array a)))


(defn make-stem-array
  ([square] (make-stem-array square square))
  ([rows cols] (make-array Double/TYPE rows cols)))

(defmacro aget!
  ([array x]      `(aget ~(vary-meta array assoc :tag 'doubles) ~x))
  ([array x & ys] `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~x)]
                     (aget! a# ~@ys))))

(defmacro aset! [array x y v]
  (let [nested-array `(aget ~(vary-meta array assoc :tag 'objects) ~x)
        a-sym         (with-meta (gensym "a") {:tag 'doubles})]
    `(let [~a-sym ~nested-array]
       (aset ~a-sym ~y (double ~v)))))
