(ns stem.util
  (:use [clojure.pprint] [stem.constants])
  (:require [clojure.string :as str]
            [clojure.java.io :as io]
            [clojure.zip :as z])
  (:import [java.io File StringReader BufferedReader FileNotFoundException]
           [org.yaml.snakeyaml Yaml]
           [java.util Random]))

(def in-production? false)

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

(defn make-tree-zipper
  "Makes a valid zipper based on the nested tree structure used throughout STEM:
  [{node info} [left branch][right branch]]"
  [tree]
  (z/zipper second rest #(vec (cons (first %1) %2)) tree))


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

(defn seq-locs-by-id
  "Returns seq of locs in tree that match up with n-ids"
  [root n-ids]
  (loop [loc (z/next (make-tree-zipper root))
         locs (transient [])]
    (if (z/end? loc)
      (persistent! locs)
      (let [[{name :name} _ _] (z/node loc)]
        (if (contains? n-ids name)
          (recur (z/next loc) (conj! locs loc))
          (recur (z/next loc) locs))))))

(defn outgroup?
  "Returns true if loc is an outgroup, i.e., root is its parent."
  [loc]
  (let [[{name :name} _ _] (z/node (z/up loc))]
    (= name *root-name*)))

(defn sisters?
  "Returns true if any of the poss-sisters are actually a sister node
  of loc."
  [loc maybe-sisters]
  (let [[{sis-name :name} _ _] (z/node (or (z/left loc) (z/right loc)))]
    (contains? maybe-sisters sis-name)))

(defn check-hybrid-loc
  [loc h-specs]
  (let [err-str "Not a valid hybrid tree. "
        [{name :name} _ _] (z/node loc)]
    (cond
     (outgroup? loc) (abort (str err-str name " is an outgroup.\n"))
     (sisters? loc h-specs) (abort (str err-str name " is a sister of another hybrid species.\n")))))

(defn is-valid-hybrid-tree
  "Checks that hybrid species don't belong to outgroup and, if there are more
  than two species, that none of them are sisters."
  [tree h-specs]
  (let [locs (seq-locs-by-id tree h-specs)]
    (doseq [l locs] (check-hybrid-loc l h-specs))))


(defn zero->min [val]
  (if-not (zero? val) val (Double/MIN_VALUE)))

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
  (format "%.5f" (double time))
  ;(cl-format nil "~,5f" time)
  )

(defn make-precise [num]
  (-> num (format-time) (to-double)))

(defn remove-whitespace [s]
  (str/replace s #"\s" ""))

(defn indexed-map
  [coll]
  (zipmap coll (range (count coll))))

;;;;;;;;;;; functions for parsing yaml settings file ;;;;;;;;;;;;;;;;;;;;;;;
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

(defn file->yaml-map
  "Build yaml map from yaml file.  Replaces all tabs with 2 spaces to ensure
  yaml parser doesn't choke."
  ([f] (file->yaml-map f ""))
  ([f e-str]
     (let [yaml (Yaml.)]
       (with-exc (.load yaml (str/replace (read-file f) "\t" "  ")) e-str))))

(defn yaml-map->map
  "Turns a yaml 'LinkedHashMap' into a clojure map, where the keys are
  keywords"
  [y-map]
  (zipmap (map keyword (.keySet y-map)) (map #(into {} %) (.values y-map))))

(defn settings-file->map
  [f e-str]
  (-> (file->yaml-map f e-str) (yaml-map->map)))

(defn parse-settings-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files = a map - keys are filenames, values are the rate constant
  :species = a map - keys are species names, values are comma separated lineages
  :properties = a map with user configurable properties.

  In addition, this function attempts to replace all tabs with 2 spaces in the
  settings file, so the yaml parser doesn't choke."
  ([] (parse-settings-file (get-settings-filename)))
  ([f] (parse-settings-file f ""))
  ([f e-str]
     (let [f-name (if f f (get-settings-filename) )
           s-map (settings-file->map f-name e-str)]
       s-map)))

(defn check-settings-map
  [m e-strs]
  (doseq [k [:properties :species]]
    (abort-if-empty (k m) (k e-strs)))
  m)

(defn parse-tree-file [f]
  (-> (with-exc (read-file f)
        (str "You must specify a species tree in '" f "' to compute the likelihood."))
      (remove-whitespace)
      (str/split #";")))

(defn get-file
  [f-name]
  (try
   (File. f-name)
   (catch Exception _
     (abort (str "Couldn't find the file: " f-name)))))

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
     (let [lid (spec-to-idx l-spec)
           rid (spec-to-idx r-spec)]
       (when (some nil? [lid rid])
         (abort (str "Either " l-spec " or " r-spec " does not match the lineages and species entered in the settings file.")))
       (get-from-upper-triangular spec-mat lid rid)))))
