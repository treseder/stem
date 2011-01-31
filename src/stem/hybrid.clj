(ns stem.hybrid
  (:use [stem.constants] [clojure.set])
  (:require [stem.newick :as n]
            [stem.lik :as l]
            [stem.util :as u]
            [clojure.zip :as z]
            [clojure.contrib.combinatorics :as c]))


(defn seq-gamma-vals
  "Returns a lazy sequence of lists, each containing possible
  gamma values for n gammas, e.g. for two gammas:
  (0 0) (0 0.01) (0 0.02)...(0.01 0) (0.01 0.01)...(1.0 1.0)"
  [n]
  (apply c/cartesian-product (repeatedly n #(range 0 1 0.01))))

(defn seq-gamma-coefs
  "Returns the seq of multiplied gammas that are the coefficients used in computing
  the likelihood for each iteration."
  [gs]
  (let [gs-with-comp (interleave gs (map #(- 1 %) gs)) ; adds 1 - gamma for each gamma
        gammas-cp (apply c/cartesian-product (partition 2 gs-with-comp))]
    (map #(apply * %) gammas-cp)))

(defn seq-liks-for-gtree
  "Returns a seq of likelihoods of a gene-tree given a seq of species trees"
  [spec-trees memo-lik-fn spec-to-lin theta g-tree]
  ; takes exp because the lik function returns the log likelihood
  (map #(Math/exp (memo-lik-fn (:vec-tree g-tree) % spec-to-lin (/ 2 theta))) spec-trees))

(defn lik-for-gammas
  "Given a seq of possible gamma values, compute the likelihood"
  [gs g-trees lik-for-gtree-fn]
  (let [g-coefs (seq-gamma-coefs gs)]
    (reduce
     (fn [total g-tree]
       (+ total
          (Math/log (->> (map * (lik-for-gtree-fn g-tree) g-coefs)
                         (reduce + 0)))))
     0
     g-trees)))

(defn find-gammas
  "Search for and return gammas that produce the largest likelihood"
  [g-trees spec-trees spec-to-lin theta num-hyb-events]
  (let [seq-gamma-vals (seq-gamma-vals num-hyb-events)
        memo-lik-fn (memoize l/calc-lik-for-tree)
        lik-gtree-fn (partial seq-liks-for-gtree spec-trees memo-lik-fn spec-to-lin theta)]
    (reduce
     (fn [[lik curr-gammas] new-gammas]
       (let [new-lik (lik-for-gammas new-gammas g-trees lik-gtree-fn)]
         (if (> new-lik lik)
           [new-lik new-gammas]
           [lik curr-gammas])))
     [Double/NEGATIVE_INFINITY (first seq-gamma-vals)]
     seq-gamma-vals)))

(defn compute-aic
  [lik k]
  (+ (* (- 2) lik) (* 2 k)))

(defn compute-k
  [hybs specs]
  (+ hybs (- specs 1)))

;; functions to build hybrid trees based on selected hybrid
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

(defn get-specs-exclude-hybs
  [hybs l-name l-descs r-name r-descs]
  ;; if hybrid is a direct child, then use it in calculating c-time, else don't
  (let [r-specs (if-not (empty? r-descs) r-descs #{r-name})
        l-specs (if-not (empty? l-descs) l-descs #{l-name})]
    (if (or (contains? hybs l-name) (contains? hybs r-name))
      [r-specs l-specs]
      [(difference r-specs hybs) (difference l-specs hybs)])))

(defn get-min-time
  [children spec-mat spec-to-idx hybs]
  (let [[[{l-name :name l-descs :desc}] [{r-name :name r-descs :desc}]] children
       [r-specs l-specs] (get-specs-exclude-hybs hybs l-name l-descs r-name r-descs)]
    (u/min-coal-time-for-node l-specs r-specs spec-mat spec-to-idx)))

(defn change-time-for-node
  "Change node time if the optimized c-time (min) is not equal to loc's c-time,
  or if the parents c-time is < loc's c-time."
  [loc epsilon spec-mat spec-idx hybs]
  (if (z/branch? loc)
    (let [min-c-time (get-min-time (z/children loc) spec-mat spec-idx hybs)]
      (if (z/up loc)
        (fix-branch loc min-c-time epsilon)
        (fix-root loc min-c-time)))
    loc))

(defn fix-tree-times
  "For each node in the tree get the min c-time, and make sure the
  molecular clock holds"
  [tree spec-matrix spec-to-index hybs]
  (loop [loc (z/zipper second rest #(vec (cons (first %1) %2)) tree)]
    (if (z/end? loc)
      (z/root loc)
      (recur (z/next (change-time-for-node loc 1e-10 spec-matrix spec-to-index hybs))))))



(defn get-fixed-descs [children]
  (let [[[{l-name :name l-descs :desc}] [{r-name :name r-descs :desc}]] children
        r-specs (if-not (empty? r-descs) r-descs #{r-name})
        l-specs (if-not (empty? l-descs) l-descs #{l-name})]
    (union r-specs l-specs)))

(defn get-sibling-loc [loc] (or (z/left loc) (z/right loc)))

(defn change-grand-parent
  [gp-loc h h-sib p p-sib]
  (let [[{name :name c-time :c-time} _ _] p
        gp-l-child [(n/create-node name 0 c-time (get-fixed-descs [h p-sib])) h p-sib]
        [{gp-name :name gp-c-time :c-time} _ _] (z/node gp-loc)
        new-gp-branch [(n/create-node gp-name 0 gp-c-time (get-fixed-descs [h-sib gp-l-child]))
                       h-sib gp-l-child]]
    (z/replace gp-loc new-gp-branch)))

(defn change-hybrid-position
  [loc]
  (let [h (z/node loc)
        h-sib (z/node (get-sibling-loc loc))
        p-loc (z/up loc)
        p-sib  (z/node (get-sibling-loc p-loc))
        gp-loc (-> p-loc z/up)]
    (z/root (change-grand-parent gp-loc h h-sib (z/node p-loc) p-sib))))

(defn permute-hybrid-tree-for-spec
  [tree spec]
  (loop [loc (z/zipper second rest #(vec (cons (first %1) %2)) tree)]
    (cond
     (z/end? loc) (z/root loc)
     (= (:name (first (z/node loc))) spec) (change-hybrid-position loc)
     :else (recur (z/next loc)))))

(defn make-hybrid-trees-for-spec
  [trees spec]
  (reduce
   #(conj %1 %2 (permute-hybrid-tree-for-spec %2 spec))
   '() 
   trees))

(defn make-hybrid-trees-for-specs
  [tree specs]
  (reduce
   #(make-hybrid-trees-for-spec %1 %2)
   [tree]
   specs))
