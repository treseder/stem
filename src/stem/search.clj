(ns stem.search
  (:use [stem.constants] [clojure.set])
  (:require  [stem.util :as u]
             [stem.messages :as m]
             [stem.lik :as l]
             [stem.newick :as n]
             [clojure.zip :as z]))

(defn change-time-for-node [loc parent-loc epsilon]
  (let [[{name :name c-time :c-time :as node desc :desc}] (z/node loc)
        [{pc-time :c-time}] (z/node parent-loc)]
    (if (or (zero? c-time) (> pc-time c-time))
      loc
      (z/replace loc (vec (cons (n/create-node name epsilon (- pc-time epsilon) desc) (z/children loc)))))))


(defn fix-tree-times
  "After a tree has been changed, it's possible that a parent node has a smaller
  coalescent time than its children, which in reality is not possible.  This function
  fixes the tree so that all descendents have slightly smaller coalescent times than
  their parent"
  [tree]
  (let [make-node #(vec (cons (first %1) %2))]
    (loop [loc (z/zipper second rest make-node tree)]
      (if (z/end? loc)
        (z/root loc)
        (recur (z/next
                (if-let [parent (z/up loc)]
                  (change-time-for-node loc parent 0.00005)
                  loc)))))))

(defn find-target-node
  [zipped-tree target-name]
  (loop [node (z/next zipped-tree)]
    (let [name (:name (first (z/node node)))]
      (if (= name target-name)
        node
        (if-not (z/end? node)
          (recur (z/next node))
          (u/abort (str "Error: traversed tree but couldn't find " target-name ".") nil))))))

(defn get-sib-node [loc]
  "Figures out if sibling node is left or right and returns it"
  (-> (or (z/left loc) (z/right loc))
      (z/node)))

(defn change-tree
  "Swaps the sibling of loc with loc's left or right child"
  [loc left-child?]
  (let [to-sib-loc (if (z/left loc) z/left z/right)
        sib-node (-> loc (to-sib-loc) (z/node))
        children (z/children loc)
        child-node (if left-child? (first children) (second children))
        to-child-loc (if left-child? z/down #(-> % (z/down) (z/right)))]
    (-> loc (to-child-loc) (z/replace sib-node)
        (z/up) (to-sib-loc) (z/replace child-node) (z/root) (fix-tree-times))))

(defn min-coal-time-for-node
  [l-specs r-specs spec-mat spec-to-idx]
  (reduce min
   (for [l-spec l-specs r-spec r-specs]
     (let [i (spec-to-idx l-spec)
           j (spec-to-idx r-spec)]
       (if (< i j)
         (u/aget! spec-mat i j)
         (u/aget! spec-mat j i))))))

(defn make-node
  "As the zipper puts the tree together, this function is called to create the new
  changed parent nodes.  Each new node needs the minimum time of the pairs of species
  from each branch, along with the new descendent set."
  [spec-mat spec-to-idx node children]
  (let [[[{l-name :name l-descs :desc lc-time :c-time}] [{r-name :name r-descs :desc rc-time :c-time}]] children
        r-specs (if-not (empty? r-descs) r-descs #{r-name})
        l-specs (if-not (empty? l-descs) l-descs #{l-name})
        c-time  (min-coal-time-for-node l-specs r-specs spec-mat spec-to-idx)
        new-node (n/create-node (:name (first node))
                                (- c-time lc-time)  
                                c-time
                                (union r-specs l-specs))]
    (vec (cons new-node children))))

(defn permute-tree
  [s-vec-tree make-node-fn num-i-nodes rand-fun]
  (let [zipped-tree (z/zipper second rest make-node-fn s-vec-tree)
        ;; chooses random internal node as target
        target-name (str *internal-node-name* "-" (+ (rand-fun :int num-i-nodes) 1))
        target-loc (find-target-node zipped-tree target-name)
        target-sib-node (get-sib-node target-loc)
        ;; if rand num = 1 changes left child, else right child
        left-child? (if (= (+ (rand-int 2) 1) 1) true false)]
    (change-tree target-loc left-child?)))

(defn keep-tree?
  "Even if the lik is less than the previous lik, still keep the tree
  based on a certain probability (to overcome local minimums)"
  [prev-lik lik i c0 beta rand-fun]
  (let [tree-prob (-> (- lik prev-lik) (/ (/ c0 (+ 1 (* i beta)))) (Math/exp))]
    (> tree-prob (rand-fun))))

(defn maybe-add-to-best
  [lik tree trees keep-n]
  (let [[min-lik min-tree :as min-entry] (last trees)]
    (cond (< (count trees) keep-n) (conj trees [lik tree])
          (> lik min-lik) (-> (conj trees [lik tree]) (disj min-entry))
          :default trees)))

(defn search-for-trees
  [s-vec-tree gene-trees spec-matrix props env]
  ;(n/see-vector-tree s-vec-tree)
  ;(u/print-array spec-matrix)
  (let [make-node-fn (partial make-node spec-matrix (:spec-to-index env))
        int-node-cnt (- (count (:spec-to-index env)) 2)
        theta (:theta env)
        beta (get props "beta" *beta-default*)
        total-iters (get props "bound_total_iter" *bound-total-iter-default*)
        burnin (get props "burnin" *burnin-default*)
        rand-fun (u/rand-generator (props "seed"))
        num-saved-trees (get props "num_saved_trees" *num-saved-trees-default*)]
   (loop [prev-lik (l/calc-mle gene-trees s-vec-tree (:spec-to-lin env) theta)
          s-tree s-vec-tree
          best-trees (-> (sorted-set-by #(> (first %1) (first %2))) (conj [prev-lik s-tree]))
          c0 (* (- prev-lik) 0.25)
          max-lik-change 0.0
          iter 1]
     (if (= iter total-iters)
       (take num-saved-trees best-trees)
       (let [tree (permute-tree s-tree make-node-fn int-node-cnt rand-fun)
             lik (l/calc-mle gene-trees tree (:spec-to-lin env) theta)
             abs-dif (Math/abs (double (- prev-lik lik)))
             next-c0 (if (> iter burnin) max-lik-change c0)]
         (if (or (> lik prev-lik) (keep-tree? prev-lik lik iter c0 beta rand-fun))
           (recur lik tree (maybe-add-to-best lik tree best-trees num-saved-trees) next-c0 (max abs-dif max-lik-change) (inc iter))
           (recur prev-lik s-tree best-trees next-c0 (max abs-dif max-lik-change) (inc iter))))))))

