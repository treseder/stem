(ns stem.search
  (:use [stem.constants] [clojure.set])
  (:require  [stem.util :as u]
             [stem.messages :as m]
             [stem.lik :as l]
             [stem.newick :as n]
             [clojure.zip :as z]))

(defn fix-root
  [loc min-c-time]
  (let [[{name :name c-time :c-time desc :desc}] (z/node loc)]
    (if (not= min-c-time c-time)
      (z/replace loc (vec (cons (n/create-node name 0 min-c-time desc) (z/children loc))))
      loc)))

(defn fix-branch
  [loc min-c-time epsi]
  (let [[{name :name c-time :c-time desc :desc}] (z/node loc)
        [{pc-time :c-time}] (z/node (z/up loc))]
    (cond
     (> min-c-time pc-time) (z/replace loc (vec (cons (n/create-node name 0 (- pc-time epsi) desc) (z/children loc))))
     (not= min-c-time c-time) (z/replace loc (vec (cons (n/create-node name 0 min-c-time desc) (z/children loc))))
     :else loc)))

(defn get-min-time
  [children spec-mat spec-to-idx]
  (let [[[{l-name :name l-descs :desc}] [{r-name :name r-descs :desc}]] children
        r-specs (if-not (empty? r-descs) r-descs #{r-name})
        l-specs (if-not (empty? l-descs) l-descs #{l-name})]
    (u/min-coal-time-for-node l-specs r-specs spec-mat spec-to-idx)))

(defn change-time-for-node
  "Change node time if the optimized c-time (min) is not equal to loc's c-time,
  or if the parents c-time is < loc's c-time."
  [loc epsilon spec-mat spec-idx]
  (if (z/branch? loc)
    (let [min-c-time (get-min-time (z/children loc) spec-mat spec-idx)]
      (if (z/up loc)
        (fix-branch loc min-c-time epsilon)
        (fix-root loc min-c-time)))
    loc))

(defn fix-tree-times
  "For each node in the tree get the min c-time, and make sure the
  molecular clock holds"
  [tree spec-matrix spec-to-index]
  (loop [loc (z/zipper second rest #(vec (cons (first %1) %2)) tree)]
    (if (z/end? loc)
      (with-meta (z/root loc) (meta tree))
      (recur (z/next (change-time-for-node loc 1e-10 spec-matrix spec-to-index))))))

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
        (z/up) (to-sib-loc) (z/replace child-node) (z/root))))

(defn make-node
  "As the zipper puts the tree together, this function is called to create the new
  changed parent nodes.  Each new node needs the minimum time of the pairs of species
  from each branch, along with the new descendent set."
  [node children]
  (let [[[{l-name :name l-descs :desc}] [{r-name :name r-descs :desc}]] children
        r-specs (if-not (empty? r-descs) r-descs #{r-name})
        l-specs (if-not (empty? l-descs) l-descs #{l-name})
        ;c-time  (u/min-coal-time-for-node l-specs r-specs spec-mat spec-to-idx)
        new-node (n/create-node (:name (first node))
                                0 ;; sets zero for now - time will be
                                0 ;; fixed later in the fix-times-fn
                                (union r-specs l-specs))]
    (vec (cons new-node children))))

(defn permute-tree
  [s-vec-tree make-node-fn num-i-nodes rand-fun spec-matrix spec-to-index]
  (let [zipped-tree (z/zipper second rest make-node-fn s-vec-tree)
        ;; chooses random internal node as target
        target-name (str *internal-node-name* "-" (+ (rand-fun :int num-i-nodes) 1))
        target-loc (find-target-node zipped-tree target-name)
        target-sib-node (get-sib-node target-loc)
        ;; if rand num = 1 changes left child, else right child
        left-child? (if (= (+ (rand-fun :int 2) 1) 1) true false)]
    (fix-tree-times (change-tree target-loc left-child?) spec-matrix spec-to-index)))

(defn keep-tree?
  "Even if the lik is less than the previous lik, still keep the tree
  based on a certain probability (to overcome local minimums)"
  [prev-lik lik i c0 beta rand-fun]
  (let [tree-prob (-> (- lik prev-lik) (/ (/ c0 (+ 1 (* i beta)))) (Math/exp))]
    (or (> lik prev-lik) (> tree-prob (rand-fun)))))

(defn maybe-add-to-best
  [lik tree trees keep-n]
  (let [min-sr (last trees)
        new-entry [lik tree]]
    (cond (< (count trees) keep-n) (conj trees new-entry)
          (>= lik (first min-sr)) (-> (conj trees new-entry) (disj min-sr))
          :default trees)))

(defn print-perc-complete
  "prints % complete to out so users have some idea
  when the search will finish"
  [percent]
  (do (print "\r") (print (str percent "%") "completed") (flush)))

(defn sort-fn
  [[lik tree] [lik2 tree2]]
  ;(swank.core/break)
  (cond
   (> lik lik2) -1
   (< lik lik2) 1
   (and (= lik lik2) (u/quasi-isomorphic? tree tree2)) 0
   :default (compare (n/vector-tree->newick-str tree true)
                     (n/vector-tree->newick-str tree2 true))))

(defn search-for-trees
  [s-vec-tree gene-trees spec-matrix props env]
  (let [int-node-cnt (- (count (:spec-to-index env)) 2)
        theta (:theta env)
        beta (get props "beta" *beta-default*)
        total-iters (get props "bound_total_iter" *bound-total-iter-default*)
        burnin (get props "burnin" *burnin-default*)
        rand-fun (u/rand-generator (props "seed"))
        num-saved-trees (get props "num_saved_trees" *num-saved-trees-default*)]
   (loop [prev-lik (l/calc-lik gene-trees s-vec-tree (:spec-to-lin env) theta)
          curr-tree s-vec-tree
          ;; sorted, descending, by likelihood
          best-trees (sorted-set-by sort-fn [prev-lik curr-tree])
          c0 (* (- prev-lik) 0.25)
          max-lik-change 0.0
          iter 1]
     (print-perc-complete (int (* (/ iter total-iters) 100)))
     (if (= iter total-iters)
       ; done -  return best trees
       (take num-saved-trees best-trees)
       ; continue permuting trees
       (let [new-tree (permute-tree curr-tree make-node int-node-cnt rand-fun spec-matrix (env :spec-to-index))
             lik (l/calc-lik gene-trees new-tree (:spec-to-lin env) theta)
             abs-dif (Math/abs (double (- prev-lik lik)))
             next-c0 (if (> iter burnin) max-lik-change c0)]
         (if (keep-tree? prev-lik lik iter c0 beta rand-fun)
           (recur lik new-tree
                  (maybe-add-to-best lik new-tree best-trees num-saved-trees) next-c0
                  (max abs-dif max-lik-change) (inc iter))
           (recur prev-lik curr-tree best-trees next-c0
                  (max abs-dif max-lik-change) (inc iter))))))))

