(ns stem.tree
  (:use [stem.constants])
  (:require [stem.util :as util]
            [clojure.set :as c-set]
            [clojure.string :as str]
            [stem.print-tree :as pt]))

(defprotocol StemTree
  (tree->seq [tree] [tree keep?])
  (tree->vec-tree [tree])
  (tree->newick [tree])
  (print-tree [tree]))

(defprotocol StemNode
  (is-leaf? [node])
  (children [node])
  (node->seq [node])
  (node->vec [node]))

(defrecord Node [^String name ^double b-len ^double c-time
                 ^clojure.lang.PersistentHashSet desc l r]
  StemNode
  (is-leaf? [node] (not (or l r))) ;; a leaf it both l and r are nil

  (node->seq
   [node]
   (if (is-leaf? node)
     node
     (lazy-cat (node->seq l) (node->seq r) node)))

  (children [node] (if (is-leaf? node) '() '(l r)))
  
  (node->vec
   [node]
   (if (is-leaf? node)
     [(str name)]
     [(str name ":" (format "%1.5f" c-time)) (node->vec l) (node->vec r)])))

(defrecord Tree [root rate gamma-vals]
  StemTree
  (tree->seq [tree keep?] (filter keep? (node->seq root)))
  
  (tree->newick [tree])
  
  (tree->vec-tree
   [tree]
   (node->vec root))
  
  (print-tree
   [tree]
   (-> (tree->vec-tree tree) (pt/draw-binary-tree))))


(defn make-node
  [name b-len c-time desc l r]
  (Node. name b-len c-time desc l r))

(defn make-tree [root rate g-vals] (Tree. root rate g-vals))


(defn max-c-time
  "Normal way to compute the coalescent time of a tree given branch lengths
  in the Newick format.
  Given a nodes left and right children nodes, calculate its c-time"  
  [l r]
  (let [{l-len :b-len l-time :c-time} l
        {r-len :b-len r-time :c-time} r]
    (max (+ l-len l-time) (+ r-len r-time))))

(defn optimized-c-time
  "Computes the coalescent times for optimized trees using the D-AB matrix, i.e., the
  min-species-matrix times"
  [spec-to-index spec-mat l r]
  (let [{l-name :name l-descs :desc} l
        {r-name :name r-descs :desc} r
        l-specs (if (is-leaf? l) #{l-name} l-descs)
        r-specs (if (is-leaf? r) #{r-name} r-descs)]
    (util/min-coal-time-for-node l-specs r-specs spec-mat spec-to-index)))

(defn merge-desc
  "Builds the descendent set of a parent given its direct descendent nodes"  
  [l r]
  (let [{l-name :name l-desc :desc} l
        {r-name :name r-desc :desc} r
         l-d (if (is-leaf? l) #{l-name} l-desc)
         r-d (if (is-leaf? r) #{r-name} r-desc)]
    (c-set/union r-d l-d)))

;; counter to uniquely identify internal nodes
(def i-node-counter (util/make-counter 0))

(defn rec-newick->tree
  "Builds the nested vector tree structure that will be used throughout
  the rest of the program.
  f is a function that takes two children nodes and determines the
  c-time of the parent."  
  [el div f]
  (if (= (count el) 2)
    ;; leafs are (name time)
    (let [[n t] el] 
      (make-node (str n) (/ t div) 0.0 #{} nil nil))
    (let [[l-el r-el t] el
          l (rec-newick->tree l-el div f)
          r (rec-newick->tree r-el div f)]
      ;; internal node
      (make-node (str *internal-node-name* "-" ((i-node-counter :next)))
                 (if-not (zero? t) (/ t div) t)
                 (f l r)
                 (merge-desc l r)
                 l r))))

(defn make-branch-lens-parseable
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
  [s]
  (str "("
       (-> s
           (str/replace ";" "")
           (str/replace "," ")(")
           (make-branch-lens-parseable))
       ")"))

(defn newick->tree
  "Parses s and builds a binary tree structure as a vector of
  vectors.  s must conform to proper Newick format.  If the
  root is named, it is ignored.  Leafs are assumed to contain a
  name:branch-len.  All interior nodes must have a branch-len,
  but not named."
  ([s] (newick->tree s 1 1 max-c-time))
  ([s rate theta] (newick->tree s rate theta max-c-time))
  ([s rate theta f]
     ((i-node-counter :reset))   ; resets counter for each tree parsed
     (let [prepped-str (make-str-parseable s)
             divisor (* rate theta)
             lst (read-string prepped-str)
             left (rec-newick->tree (first lst) divisor f)
             right (rec-newick->tree (second lst) divisor f)
             c-time (f left right)
             desc-set (merge-desc left right)]
         (make-tree (make-node *root-name* 0.0 c-time desc-set left right) rate nil))))

