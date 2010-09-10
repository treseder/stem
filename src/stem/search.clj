(ns stem.search
  (:use [stem.constants] [clojure.set])
  (:require  [stem.util :as u]
             [stem.messages :as m]
             [stem.lik :as l]
             [stem.newick :as n]
             [clojure.zip :as z]))


(defn find-target-node
  [zipped-tree target-name]
  (println (str "looking for target: " target-name))
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
        (z/up) (to-sib-loc) (z/replace child-node) (z/root))))

(defn min-coal-time-for-node
  [l-specs r-specs spec-mat spec-to-idx]
  (reduce min
   (for [l-spec l-specs r-spec r-specs]
     (let [i (l-spec spec-to-idx)
           j (r-spec spec-to-idx)]
       (if (< i j)
         (u/aget! spec-mat i j)
         (u/aget! spec-mat j i))))))

(defn make-node
  "As the zipper puts the tree together, this function is called to create the new
  changed parent nodes.  Each new node needs the minimum time of the pairs of species
  from each branch, along with the new descendent set."
  [spec-mat spec-to-idx node children]
  (println (str "making node with... " node))
  (let [r-specs (:desc (first (second children)))
        l-specs (:desc (first (first children)))
        new-node (n/create-node (:name (first node))
                                (min-coal-time-for-node l-specs r-specs spec-mat spec-to-idx)
                                0.0
                                (union r-descs l-descs))]
      (vec (cons new-node children))))

(defn permute-tree
  [s-vec-tree make-node-fn num-i-nodes]
  (let [zipped-tree (z/zipper second rest make-node-fn s-vec-tree)
        ;; chooses random internal node as target
        target-name (str *internal-node-name* "-" (+ (rand-int num-i-nodes) 1))
        target-loc (find-target-node zipped-tree target-name)
        target-sib-node (get-sib-node target-loc)
        ;; if rand num = 1 changes left child, else right child
        left-child? (if (= (+ (rand-int 2) 1) 1) true false)]
    (println (str "Left child: " left-child?))
    (change-tree target-loc left-child?)))

(defn keep-tree?
  "Even if the lik is less than the previous lik, still keep the tree
  based on a certain probability to get out of local minimums"
  [prev-lik lik i c0 beta]
  (-> (- lik prev-lik)
      (/ (/ c0 (+ 1 (* i beta))))
      (Math/exp)
      (> (rand))))

(defn maybe-add-to-best
  [lik tree trees keep-n]
  (let [[min-lik min-tree :as min-entry] (first (last trees))]
    (cond (< (count trees) keep-n) (conj! trees [lik tree])
          (> lik min-lik) (-> (conj! trees [lik tree]) (disj! min-entry))
          trees)))

(defn search-for-trees
  [s-vec-tree gene-trees spec-to-lin spec-matrix spec-to-idx props iters num-trees-to-save]
  (let [make-node-fn (partial make-node spec-matrix spec-to-idx)
        int-node-cnt (- (count spec-to-idx) 2)
        theta (:theta props)
        beta (:beta props)]
   (loop [prev-lik (l/calc-mle gene-trees s-vec-tree spec-to-lin theta)
          s-tree s-vec-tree
          best-trees (transient (sorted-set))
          c0 (* (- prev-lik) 0.25)
          max-lik-change
          iter 1]
     (if (= iter iters)
       (persistent! best-trees)
       (let [tree (permute-tree s-tree make-node-fn int-node-cnt)
             lik (l/calc-mle gene-trees tree spec-to-lin theta)
             abs-dif (Math/abs (- prev-lik lik))
             next-c0 (if (> iter 200) max-lik-change c0)]
         (if (or (> lik prev-lik) (keep-tree? prev-lik lik iter c0 beta))
           (recur lik tree (maybe-add-to-best lik tree best-trees) next-c0 (max abs-dif max-lik-change) (inc iters))
           (recur prev-lik s-tree best-trees next-c0 (max abs-dif max-lik-change) (inc iters))))))))

