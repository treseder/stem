(ns stem.search
  (:use [stem.constants])
  (:require  [stem.util :as u]
             [stem.messages :as m]
             [stem.lik :as lik]
             [clojure.zip :as z]))


(defn find-target-node
  [zipped-tree target-name]
  (loop [node (z/next zipped-tree)]
    (let [name (:name (first (z/node node)))]
      (if (= name target-name)
        node
        (if-not (z/end? node)
          (recur (z/next node))
          (u/abort (str "Error: traversed tree but couldn't find " name ".") nil))))))

(defn get-sib-node [loc]
  "Figures out if sibling node is left or right and returns it"
  (-> (or (z/left loc) (z/right loc))
      (z/node)))

(defn change-tree
  "Swaps the sibling of loc with locs left or right child"
  [loc left-child?]
  (let [to-sib-loc (if (z/left loc) z/left z/right)
        sib-node (-> loc (to-sib-loc) (z/node))
        children (z/children loc)
        child-node (if left-child? (first children) (second children))
        to-child-loc (if left-child? z/down #(-> % (z/down) (z/right)))]
    (-> loc (to-child-loc) (z/replace sib-node)
        (z/up) (to-sib-loc) (z/replace child-node) (z/root))))

(defn make-node [node children]
  (vec (cons (first node) children)))

(defn permute-tree
  [s-vec-tree num-i-nodes]
  (let [zipped-tree (z/zipper second rest make-node s-vec-tree)
        ;; chooses random internal node as target
        target-name (str *internal-node-name* "-" (+ (rand-int num-i-nodes) 1))
        target-loc (find-target-node zipped-tree target-name)
        target-sib-node (get-sib-node target-loc)]
    (case (+ (rand-int 3) 1)
          ;; changes sibling and left child
          1 (change-tree target-loc true)
          ;; changes sibling and right child
          2 (change-tree target-loc false)
          ;; no changes tree remains the same
          3 (z/root s-vec-tree))))

(defn search-for-trees
  [s-vec-tree gene-trees num-lins iters]
  
  )

