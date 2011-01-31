(ns stem.util
  (:use [clojure.pprint] [stem.constants])
  (:require [clojure.string :as str]
            [clojure.java.io :as io])
  (:import [java.io File StringReader BufferedReader FileNotFoundException]
           [java.util Random]))


(defn quasi-isomorphic?
  "Returns true if the two trees are quasi-isomorphic *and* the contents
  of the leaves match"
  [t1 t2]
  (let [[{name1 :name} l1 r1] t1
        [{name2 :name} l2 r2] t2]
    (cond
     ; comparing two leaves
     (every? nil? [l1 l2 r1 r2]) (= name1 name2)
     ; compares a branch and leaf
     (some nil? [l1 l2 r1 r2]) false
     ; compares two branches: either the left nodes of each are equal
     ; or the left and right of each are equal
     :default (or
               (and (quasi-isomorphic? l1 l2) (quasi-isomorphic? r1 r2))
               (and (quasi-isomorphic? l1 r2) (quasi-isomorphic? r1 l2))))))

(def in-production? false)

(defmacro with-exc [body message]
  `(try
    ~body
    (catch Exception e# (stem.util/abort ~message e#))))

(defn rand-generator
  "Returns a random number generator given a seed."
  [& [seed]]
  (let [rg (if seed (Random. seed) (Random.))]
    (fn [& [type n]] (if (= type :int) (.nextInt rg n) (.nextDouble rg)))))

(defn zero->tiny-num
  "Zero times aren't really valid, but sometimes users include them; change
  to very small number "
  [num]
  (if-not (zero? num) num 0.00001))

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

(defn abort
  ([message] (abort message nil))
  ([message e]
     (println message)
     (when e
       (println (.getMessage e)))
     (when in-production?
       (do
         (println "\nExiting STEM...")
         (System/exit 1)))))

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

(defmethod to-double java.lang.Integer [i] (double i))

(defmethod to-double java.lang.Double [i] i)

(defn read-file [f]
  (try
    (slurp f)
    (catch FileNotFoundException e
      (abort (str "File " f " not found.")))))

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

(defn make-precise [num]
  (-> num (format-time) (to-double)))

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

(defn get-from-upper-triangular
  "Assumes matrix is upper triangular"
  [matrix i j]
  (if (< i j)
    (aget! matrix i j)
    (aget! matrix j i)))

(defn min-coal-time-for-node
  [l-specs r-specs spec-mat spec-to-idx]
  (reduce
   min
   (for [l-spec l-specs r-spec r-specs]
     (get-from-upper-triangular spec-mat (spec-to-idx l-spec) (spec-to-idx r-spec)))))
