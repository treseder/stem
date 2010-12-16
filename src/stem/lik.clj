(ns stem.lik
  (:use [clojure.contrib pprint])
  (:require [stem.newick :as newick]
            [clojure.set :as set]
            [stem.gene-tree :as g-tree]
            [stem.util :as util]
            [clojure.contrib.combinatorics :as comb]))


(defn calc-lik-for-coalescent-event
  [num-lins two-div-theta start end leaving?]
  (let [comb (* num-lins (- num-lins 1))
        time-dif (util/zero->tiny-num (- end start))
        exp-part (Math/exp (- (* comb time-dif)))
        mle  (if-not leaving? (* 2 exp-part) exp-part)
        ret-mle (Math/log (if (zero? mle) (Double/MIN_VALUE) mle))]
    ret-mle))

(defn calc-lik-for-branch
  "TODO: what about c-events that occur at the exact same time?"
  [coal-nodes lins-in-branch num-lins-at-start start end two-div-theta]
  ;; the coal-filter-fn could be refactored out to avoid duplication,
  ;; but who can resist that nice closure
  (let [coal-filter-fn (fn [{c-time :c-time c-lins :desc}]
                         (and (>= c-time start)
                              (< c-time end)
                              (set/subset? c-lins lins-in-branch)))
        ;; events that occurred within branch
        coal-events (filter coal-filter-fn coal-nodes) 
        lik-fn (fn [[mle num-lins last-c-time] event]
                 [(+ mle (calc-lik-for-coalescent-event num-lins two-div-theta last-c-time (:c-time event) false)),
                  (- num-lins 1),
                  (:c-time event)])
        [cum-mle rem-num-lins last-c-time] (reduce lik-fn [0 num-lins-at-start start] coal-events)]
    ;; at the end of each branch, the probability that the remaining
    ;; lineages don't coalesce must be calculated
    (let [ret-mle (if (> rem-num-lins 1)
                    [(+ cum-mle
                        (calc-lik-for-coalescent-event
                         rem-num-lins two-div-theta last-c-time end true)), rem-num-lins]
                    [cum-mle, rem-num-lins])]
      ret-mle)))

(defn calc-lik-for-branches
  "Calculates mle for each branch and returns a vector of the form:
  [cum-mle free-lins start-time].  Each branch mle is computed in the
  parent node because the calc needs the stop time of the molecular
  clock for that particular branch"
  [s-node coal-nodes spec-to-lin two-div-theta]
  (let [[n l r] s-node]
    (if-not l ;leaf
      (let [cum-lins (spec-to-lin (:name n))]
        [0 cum-lins (count cum-lins) 0])
      (let [end-time (:c-time n)
            ;; depth traversal - gets to bottom of tree and works up
            [l-cum-mle l-lins l-num-lins l-time] (calc-lik-for-branches l coal-nodes spec-to-lin two-div-theta)
            [r-cum-mle r-lins r-num-lins r-time] (calc-lik-for-branches r coal-nodes spec-to-lin two-div-theta)
            ;; calculating the two branches from this node to its two descendents
            [l-mle l-end-lins] (calc-lik-for-branch coal-nodes l-lins l-num-lins l-time end-time two-div-theta)
            [r-mle r-end-lins] (calc-lik-for-branch coal-nodes r-lins r-num-lins r-time end-time two-div-theta)
            cum-mle (+ l-cum-mle r-cum-mle l-mle r-mle)]
        [cum-mle (set/union l-lins r-lins) (+ l-end-lins r-end-lins) end-time]))))


(defn tree->sorted-internal-nodes
  "All nodes of the tree except the leaves are returned, ordered by coalescent time.
  These nodes represent the coalescent events of the gene tree."
  [tree]
  (sort-by #(:c-time %)
           (filter #(not-empty (:desc %)) (newick/tree->seq tree))))

(def counter (let [count (ref 0)] #(dosync (alter count inc))))

(defn calc-lik-for-tree
  "Calculates the maximum likelihood estimate given a
  gene tree matrix and a species tree"
  [gene-tree s-tree spec-to-lin two-div-theta]
  (let [coal-nodes (tree->sorted-internal-nodes gene-tree)
        [tree-mle coaled-lins num-lins end-time] (calc-lik-for-branches s-tree coal-nodes spec-to-lin two-div-theta)
        ;; there could still be coalescent events that occur *before* (or
        ;; after, depending on outlook) the first speciation event
        [before-tree-mle _] (calc-lik-for-branch coal-nodes coaled-lins
                                                 num-lins end-time
                                                 Double/POSITIVE_INFINITY
                                                 two-div-theta)]
;    (println (str "tree lik: " (+ tree-mle before-tree-mle)))
    (+ tree-mle before-tree-mle)))

(defn calc-lik
  [gene-trees s-tree spec-to-lin theta]
  (reduce #(+ %1 (calc-lik-for-tree (:vec-tree %2) s-tree spec-to-lin (/ 2 theta))) 0 gene-trees))

(defn calc-liks
  [gene-trees s-trees spec-to-lin theta]
  (map #(calc-lik gene-trees % spec-to-lin theta) s-trees))
