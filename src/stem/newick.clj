(ns stem.newick
  (:require [clojure.contrib.str-utils2 :as s]
            [stem.util :as util]
            [vijual :as vij]))

(defstruct node  :name :time :progeny-map)

(def *internal-node-name* "internal")

(defn create-node [n time p-map]
  (struct node n time p-map))

(defn build-tree [el]
  (if (= (count el) 1)
    (let [[n t-str] (s/split (name (first el)) #":")]
      [(create-node n (util/to-float t-str) {})])
    [(create-node *internal-node-name* (util/to-float (nth el 2)) {})
     (build-tree (nth el 1))
     (build-tree (nth el 0))]))

(defn prep-newick-str
  "To make the tree easier for parsing the commas are replaced with ')('."
  [n-str]
  (str "(" (s/replace n-str "," ")(") ")"))

(defn build-tree-from-newick-str
  "Parses s and builds a binary tree structure as a vector of
  vectors.  s must conform to proper newick format.  If the
  root is named, it is ignored.  Leafs are assumed to contain a
  name:number, where number is optional.  All interior nodes must
  have a number, but need not be named."
  [s]
  (try
   (let [prepped-str (prep-newick-str s)
         lst (read-string prepped-str)]
     [(create-node "g-root" 0 {}) (build-tree (first lst)) (build-tree (second lst))])
   (catch Exception e
     (util/abort "An error occured parsing the newick string" e))))


;;;;;;;;;;;;;;;;;;
;; drawing trees
;;;;;;;;;;;;;;;;;;

(defn build-name-from-map [m]
  (str (:name m) ":" (:time m)))

(defn create-drawable-tree [node]
  (let [[m left right] node
        d-name (build-name-from-map m)]
    (if (and (nil? left) (nil? right))
      [d-name]
      [d-name (create-drawable-tree left) (create-drawable-tree right)])))

(defn print-tree [vec-tree]
  (vij/draw-binary-tree (create-drawable-tree vec-tree)))

(def t-str "((one:1,two:2):0.2,five:5)")
