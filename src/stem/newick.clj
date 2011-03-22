(ns stem.newick
  "Provides functions that parse newick formatted strings and build a binary
  tree with the following structure:
  [{:name :b-len :c-time :desc} [... left ...] [... right ...]]

  :name name of node (leaf names from newick, internal nodes are 'internal')
  :b-len branch length - float of elapsed time to ancestor (obtained from newick)
  :descendents a set of all node's descendents

  The library also contains a function to draw an ascii representation of the
  tree"
  (:use     [clojure.pprint] [stem.constants])
  (:require [clojure.string :as str]
            [stem.util :as util]
            [clojure.set :as c-set]
            [vijual :as vij]))

;; :name is a string, :time and :c-time are floats, and :desc is a set
(defrecord Node [^String name ^double b-len ^double c-time ^clojure.lang.PersistentHashSet desc])

;; counter to uniquely identify internal nodes
(def i-node-counter (util/make-counter 0))

(defn tree->seq  "Takes a tree and returns a seq of all the nodes in depth-first order."
  [[n l r]]
  (if-not l
    [n]
    (lazy-cat (tree->seq l) (tree->seq r) [n])))

(defn create-node [n b-len c-time desc-set]
  (Node. (.trim ^String n) (util/make-precise b-len) (util/make-precise c-time) desc-set))

(defn is-leaf? [node]
  (let [[n l r] node] (nil? l)))

(defn max-c-time
  "Normal way to compute the coalescent time of a tree given branch lengths
  in the Newick format.
  Given a parent's left and right children nodes, calculate its c-time"  
  [n1 n2]
  (let [[{l-len :b-len l-time :c-time} _ _] n1
        [{r-len :b-len r-time :c-time} _ _] n2]
    (max (+ l-len l-time) (+ r-len r-time))))

(defn optimized-c-time
  "Computes the coalescent times for optimized trees using the D-AB matrix, i.e., the
  min-species-matrix times"
  [spec-to-index spec-mat left right]
  (let [[{l-name :name l-descs :desc} _ _] left
        [{r-name :name r-descs :desc} _ _] right
        l-specs (if (empty? l-descs) #{l-name} l-descs)
        r-specs (if (empty? r-descs) #{r-name} r-descs)]
    (util/min-coal-time-for-node l-specs r-specs spec-mat spec-to-index)))

(defn merge-desc
  "Builds the descendent set of a parent given its direct descendent nodes"  
  [n1 n2]
  (let [[node1 left1 _] n1 ;; add children of node
        [node2 left2 _] n2
        c1 (if (nil? left1) #{(:name node1)} #{})
        c2 (if (nil? left2) #{(:name node2)} #{})]
    (c-set/union (:desc node1) (:desc node2) c1 c2)))

(defn build-tree
  "Builds the nested vector tree structure that will be used throughout
  the rest of the program.
  c-time-fn is a function that takes two children nodes and determines the
  c-time of the parent."  
  [el div c-time-fn]
  (if (= (count el) 2)
    (let [[n t] el] ;; leafs are name:time
      ;; leaf node
      [(create-node (str n) (/ t div) 0.0 #{})])
    (let [[l r t] el
          left (build-tree l div c-time-fn)
          right (build-tree r div c-time-fn)
          c-time (c-time-fn left right)
          desc-set (merge-desc left right)]
      ;; internal node
      [(create-node (str *internal-node-name* "-" ((i-node-counter :next)))
                    (if-not (zero? t) (/ t div) t)
                    c-time
                    desc-set)
       left right])))

(defn add-branch-lens 
  "In some cases the Newick str will not include branch lens, but the
  parser assumes it does.  Adds zero branch lens to allow code reuse.
  Assumes the commas have already been replaced by ')('."
  [s]
  ;; if string contains ":" then branch lens are present
  (if (.contains s ":")
    (str/replace s ":" " ")
    (str/replace s ")" " 0.0)")))

(defn make-str-parseable
  "Makes the Newick str parseable by read-string."
  [n-str]
  (str "("
       (-> n-str
           (str/replace ";" "")
           (str/replace "," ")(")
           (add-branch-lens))
       ")"))

(defn build-tree-from-newick-str
  "Parses s and builds a binary tree structure as a vector of
  vectors.  s must conform to proper Newick format.  If the
  root is named, it is ignored.  Leafs are assumed to contain a
  name:branch-len.  All interior nodes must have a branch-len,
  but not named."
  ([s rate theta] (build-tree-from-newick-str s rate theta max-c-time))
  ([s rate theta c-time-fn]
     ((i-node-counter :reset))   ; resets counter for each tree parsed
     (try
       (let [prepped-str (make-str-parseable s)
             divisor (* rate theta)
             lst (read-string prepped-str)
             left (build-tree (first lst) divisor c-time-fn)
             right (build-tree (second lst) divisor c-time-fn)
             c-time (c-time-fn left right)
             desc-set (merge-desc left right)]
         [(create-node *root-name* 0.0 c-time desc-set) left right])
       (catch Exception e
         (util/abort "An error occured parsing the newick string" e)))))
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions to generate newick-str from tree ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn rec-tree->newick-str [node p-time branch-str]
  (if (is-leaf? node)
    (str (first node) ":"  (util/format-time p-time) branch-str)
    (let [[time l r] node]
      (str "(" (rec-tree->newick-str l time branch-str) ","
           (rec-tree->newick-str r time branch-str) ")" ":"
           (util/format-time (- p-time time))
           branch-str))))

(defn tree->newick-str 
  "Turns the likelihood tree data structure:

  [time [time [name] [name]] [name] ...]

  into a newick tree.   branch-str is any str thats printed after
  branch lengths - mainly to output trees to be used by bootstrapping in R."
  ([tree] (tree->newick-str tree ""))
  ([tree branch-str]
  (let [[n l r] tree]
    (str "(" (rec-tree->newick-str l n branch-str) ","
         (rec-tree->newick-str r n branch-str) ");" ))))

(defn sorted-branches
  [b1 b2]
  (let [[{name1 :name} _ _] b1
        [{name2 :name} _ _] b2]
    (if (pos? (compare name1 name2))
      [b1 b2]
      [b2 b1])))

(defn rec-vector-tree->newick-str
  [branch p-time]
  (let [[{:keys [name c-time]} l r] branch]
    (if (is-leaf? branch)
      (str name ":"  (util/format-time p-time))
      (let [s-brans (sorted-branches l r)]
          (str "(" (rec-vector-tree->newick-str (first s-brans) c-time) ","
            (rec-vector-tree->newick-str (second s-brans) c-time) ")" ":"
            (util/format-time (- p-time c-time)))))))

(defn vector-tree->newick-str
  "Turns any vector nested tree into a newick string.  Two trees that are
  quasi-isomorphic (with same leaf nodes) will return the same newick string."
  [tree]
  (let [[{c-time :c-time} l r] tree
        ; order branches, so quasi-isomorphic trees are the same
        s-brans (sorted-branches l r)]
    (str "(" (rec-vector-tree->newick-str (first s-brans) c-time) ","
         (rec-vector-tree->newick-str (second s-brans) c-time) ");" )))

;;;;;;;;;;;;;;;;;;
;; drawing trees
;;;;;;;;;;;;;;;;;;

(defn node-to-str [n]
  (let [[m left right] n]
    (str (str (:name m) ":") (format "%1.5f" (:c-time m)))))

(defn create-drawable-tree [node]
  (let [[m left right] node
        d-name (node-to-str node)]
    (if (and (nil? left) (nil? right))
      [d-name]
      [d-name (create-drawable-tree left) (create-drawable-tree right)])))

(defn see-vector-tree
  [vec-tree]
  (vij/draw-binary-tree (create-drawable-tree vec-tree)))

(defn see-newick-tree [n-str]
  (-> n-str (build-tree-from-newick-str 1 1) (see-vector-tree)))

(defn see-tree [vec-tree]
  (vij/draw-binary-tree vec-tree))
