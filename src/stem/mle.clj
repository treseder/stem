(ns stem.mle
  (:use [clojure.contrib pprint] [clojure set])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as util]
            [clojure.contrib.combinatorics :as comb]))


(defn subset? 
  "Is set1 a subset of set2?"
  [set1 set2]
  (and (<= (count set1) (count set2))
       (every? set2 set1)))

(defn zero->tiny-num [num]
  (if-not (zero? num) num 0.00001))

(defn factorial [n] (reduce * (range 2 (inc n))))

(defn a-choose-b [a b]
  "Returns the number of ways of choosing b items from a items."
  (/ (factorial a)
     (* (factorial b)
        (factorial (- a b)))))

(defn calc-mle-for-coalescent-event
  [num-lins two-div-theta start end leaving?]
  (let [comb (a-choose-b num-lins 2)
        exp-part (Math/exp (- (* comb (- end start) two-div-theta)))]
    (if-not leaving? (* two-div-theta exp-part) exp-part)))

(defn calc-mle-for-branch
  "TODO: what about c-events that occur at the exact same time?"
  [coal-nodes lins-in-branch num-lins-at-start start end two-div-theta]
  ;; the coal-filter-fn could be refactored out to avoid duplication,
  ;; but who can resist that nice closure
  (let [coal-filter-fn (fn [{c-time :c-time c-lins :desc}]
                         (and (>= c-time start) (< c-time end)
                              (subset? c-lins lins-in-branch)))
        coal-events (filter coal-filter-fn coal-nodes)
        mle-fn (fn [[mle num-lins last-c-time] event]
                 [(* mle (calc-mle-for-coalescent-event num-lins two-div-theta last-c-time (:c-time event) false)),
                  (- num-lins 1),
                  (:c-time event)])
        [cum-mle rem-num-lins last-c-time] (reduce mle-fn [1 num-lins-at-start start] coal-events)]
    ;; at the end of each branch, the probability that the remaining
    ;; lineages don't coalesce must be calculated, too
    (if (> rem-num-lins 1)
      [(* cum-mle (calc-mle-for-coalescent-event rem-num-lins two-div-theta last-c-time end true)), rem-num-lins]
      [cum-mle, rem-num-lins])))

(defn calc-mle-for-branches
  "Calculates mle for each branch and returns a vector of the form:
  [cum-mle free-lins start-time].  Each branch mle is computed in the
  parent node because the calc needs the stop time of the molecular
  clock for that particular branch"
  [s-node coal-nodes spec-to-lin two-div-theta]
  (let [[n l r] s-node]
    (if-not l ;leaf
      (let [cum-lins (spec-to-lin (:name n))]
        [1 cum-lins (count cum-lins) 0])
      (let [end-time (:c-time n)
            [l-cum-mle l-lins l-num-lins l-time] (calc-mle-for-branches l coal-nodes spec-to-lin two-div-theta)
            [r-cum-mle r-lins r-num-lins r-time] (calc-mle-for-branches r coal-nodes spec-to-lin two-div-theta)
            ;; calculating the two branches from this node to its two descendents
            [l-mle l-end-lins] (calc-mle-for-branch coal-nodes l-lins l-num-lins l-time end-time)
            [r-mle r-end-lins] (calc-mle-for-branch coal-nodes r-lins r-num-lins r-time end-time)]
        [(* l-cum-mle r-cum-mle) (union l-lins r-lins) (+ l-end-lins r-end-lins) end-time]))))


(defn tree->sorted-internal-nodes
  "All nodes of the tree except the leaves are returned, ordered by coalescent time.
  These nodes represent the coalescent events of the gene tree."
  [tree]
  (sort-by #(% :c-time)
           (filter #(not-empty (% :desc)) (newick/tree->seq tree))))

(defn calc-mle-for-tree
  "Calculates the maximum likelihood estimate given a
  gene tree matrix and a species tree"
  [gene-tree s-tree spec-to-lin two-div-theta]
  (let [coal-nodes (tree->sorted-internal-nodes gene-tree)
        [tree-mle coaled-lins num-lins end-time] (calc-mle-for-branches s-tree coal-nodes spec-to-lin two-div-theta)
        ;; there could still be coalescent events that occur *before* (or
        ;; after, depending on outlook) the first speciation event
        [before-tree-mle _] (calc-mle-for-branch coal-nodes coaled-lins
                                                 num-lins end-time
                                                 Double/POSITIVE_INFINITY
                                                 two-div-theta)]
    (* tree-mle before-tree-mle)))

(defn calc-mle
  [gene-trees s-tree spec-to-lin theta]
  (reduce #(+ %1 (Math/log (calc-mle-for-tree %2 s-tree spec-to-lin (/ 2 theta)))) 0 gene-trees))
