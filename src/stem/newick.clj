(ns stem.newick
  "Provides functions that parse newick formatted strings and builds a binary
  tree with the following structure:
  [{:name :time :descendents} [...] [...]]

  :name name of node (leaf names from newick, internal nodes are 'internal')
  :time float of elapsed time to ancestor (obtained from newick)
  :descendents a set of all node's descendents

  The library also contains a function to draw an ascii representation of the
  tree"
  (:require [clojure.contrib.str-utils2 :as s]
            [stem.util :as util]
            [clojure.set :as c-set]
            [vijual :as vij]))

(def *internal-node-name* "internal")
(def *root-name* "stem-root")

;; :name is a string, :time is a float, and :descendents is a ref set
(defstruct node :name :time :c-time :desc)

(defn create-node [n time c-time desc-set]
  (struct node n time c-time desc-set))

(defn max-c-time [n1 n2]
  (let [[node1 _ _] n1
        [node2 _ _] n2
        t1 (+ (:time node1) (:c-time node1))
        t2 (+ (:time node2) (:c-time node2))]
    (max t1 t2)))

(defn merge-desc [n1 n2]
  (let [[node1 left1 _] n1 ;; add children of node
        [node2 left2 _] n2
        c1 (if (nil? left1) #{(:name node1)} #{})
        c2 (if (nil? left2) #{(:name node2)} #{})]
    (c-set/union (:desc node1) (:desc node2) c1 c2)))

(defn build-tree [el]
  (if (= (count el) 1)
    (let [[n t-str] (s/split (name (first el)) #":")] ;; leafs are name:time
      ;; leaf node
      [(create-node n (util/to-float t-str) 0 #{})])
    (let [left (build-tree (nth el 1))
          right (build-tree (nth el 0))
          c-time (max-c-time left right)
          desc-set (merge-desc left right)]
      (println (str c-time " : " desc-set))
      (flush)
      ;; internal node
      [(create-node *internal-node-name*
                    (util/to-float (nth el 2))
                    c-time
                    desc-set)
       left right])))

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
         lst (read-string prepped-str)
         left (build-tree (first lst))
         right (build-tree (second lst))
         c-time (max-c-time left right)
         desc-set (merge-desc left right)]
     [(create-node *root-name* 0 c-time desc-set) left right])
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

(defn see-tree [vec-tree]
  (vij/draw-binary-tree (create-drawable-tree vec-tree)))

(def t-str "((one:1,two:2):0.2,five:5)")
