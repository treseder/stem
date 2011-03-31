(ns stem.hybrid
  (:use [stem.constants] [clojure.set])
  (:require [stem.newick :as n]
            [stem.lik :as l]
            [stem.util :as u]
            [clojure.zip :as z]
            [clojure.contrib.combinatorics :as c]))

(defrecord HybridTreeData [tree newick g-vals g-map lik k aic])

(defn create-hybrid-tree-data
  [tree newick g-vals g-map lik k aic]
  (HybridTreeData. tree newick g-vals g-map lik k aic))

(def *min-big-decimal* (BigDecimal/valueOf (Double/MIN_VALUE)))
(def *min-log* (Math/log (Double/MIN_VALUE)))
(def *math-cxt* (java.math.MathContext. 5))

(defn exp-with-prec
  [val]
  (let [cutoff-exp (int *min-log*)
        exp-val (Math/exp (/ val (Math/abs cutoff-exp)))]
    ;using valueOf to get exact representation
    (.pow (BigDecimal/valueOf exp-val) (Math/abs cutoff-exp) *math-cxt*)))

(defn log-with-prec
  "big-dec could be too small for log() resutling in -Inf, so to make it 'bigger',
   use the property that log(a) = nlog(c) + log(a/c^n)"
  [big-dec]
  (let [n (int (+ (/ (.scale big-dec) (.scale *min-big-decimal*)) 1))
        c (Double/MIN_VALUE)
        nlogc (* (Math/log c) n)
        cn (.pow *min-big-decimal* n *math-cxt*)
        a-div-cn (.divide big-dec cn *math-cxt*)]
    (+ nlogc (Math/log a-div-cn))))

(defn seq-gamma-vals
  "Returns a lazy sequence of lists, each containing possible
  gamma values for n gammas, e.g. for two gammas:
  (0 0) (0 0.01) (0 0.02)...(0.01 0) (0.01 0.01)...(1.0 1.0)"
  [n]
  (apply c/cartesian-product (repeatedly n #(range 0 1 0.01))))

(defn seq-gamma-coefs
  "Returns the seq of multiplied gammas that are the coefficients used in computing
  the likelihood for each iteration of genetrees given the species trees."
  [gs]
  (let [gs-with-compls (interleave gs (map #(- 1 %) gs)) ; adds 1 - gamma for each gamma
        gammas-cp (apply c/cartesian-product (partition 2 gs-with-compls))]
    (map #(apply * %) gammas-cp)))

(defn seq-liks-for-gtree
  "Returns a seq of likelihoods of a gene-tree given a seq of species trees"
  [spec-trees memo-lik-fn spec-to-lin theta g-tree]
  ; takes exp because the lik function returns the log likelihood
  (map #(memo-lik-fn (:vec-tree g-tree) % spec-to-lin (/ 2 theta)) spec-trees))

(defn needs-arb-prec? [lik] (< lik *min-log*))

(defn prob-for-gamma-coefs
  [liks coefs]
  (let [sum-with-prec-fn #(.add %2 %1 *math-cxt*)]
    (if (not-any? needs-arb-prec? liks)
      ; values aren't too small, use regular math functions
      (->> (map #(* (Math/exp %1) %2) liks coefs) (reduce + 0) (Math/log))
      ; uses BigDecimals for arbitrary precision
      (->> (map #(.multiply (exp-with-prec %1) (BigDecimal/valueOf (double %2)) *math-cxt*)
                liks coefs)
           (reduce sum-with-prec-fn (BigDecimal/ZERO))
           (log-with-prec)))))

(defn lik-for-gammas
  "Given a seq of possible gamma values, find the product (over all genetrees), of the sum
  of the likelihoods of each genetree given the gamma coefficients and species trees"
  [gs g-trees liks-for-gtree-fn]
  (let [g-coefs (seq-gamma-coefs gs)]
    (reduce
     (fn [total g-tree]
       (let [liks-for-gtree (liks-for-gtree-fn g-tree)
             prob (prob-for-gamma-coefs liks-for-gtree g-coefs)]
         (+ total prob)))
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
         ;(println "Lik + gammas" new-lik new-gammas)
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
  (loop [loc (u/make-tree-zipper tree)]
    (if (z/end? loc)
      ;; make sure new tree has same meta as old
      (with-meta (z/root loc) (meta tree))
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
  "Flips the hybrid node to the other possible position in the tree."
  [loc]
  (let [h (z/node loc)
        h-sib (z/node (get-sibling-loc loc))
        p-loc (z/up loc)
        p-sib  (z/node (get-sibling-loc p-loc))
        gp-loc (-> p-loc z/up)]
    (z/root (change-grand-parent gp-loc h h-sib (z/node p-loc) p-sib))))

(defn permute-tree-at-spec
  [tree spec]
  (loop [loc (u/make-tree-zipper tree)]
    (cond
     (z/end? loc) (z/root loc)
     ;; is this a hybrid node
     (= (:name (first (z/node loc))) spec) (change-hybrid-position loc)
     :else (recur (z/next loc)))))

(defn left-or-right-child
  "Checks if node's children that match name is a left or a right node.
  Returns {name 0} for left, {name 1} for right"
  [names [ _ [{l-name :name} _ _] [{r-name :name} _ _]]]
  (cond
   (contains? names l-name) {l-name 0}
   (contains? names  r-name) {r-name 1}
   :else nil))

(defn add-gamma-meta-to-tree
  "As a convention, gamma = 0 when the hybrid is a left child, and gamma = 1 when
  the hybrid node is a right child.  This function returns a map where the hybrid name
  is the key and the value is a 1, or 0.  This meta-data will be used later on to
  compute the gamma values and likelihood of hybrid trees."
  [tree specs]
  (loop [loc (u/make-tree-zipper tree)
         gamma-meta {}]
    (if (z/end? loc)
      (with-meta tree {:gamma gamma-meta})
      (if (z/branch? loc)
        (recur (z/next loc) (merge gamma-meta (left-or-right-child specs (z/node loc)))) 
        (recur (z/next loc) gamma-meta)))))

(defn add-gamma-meta-to-trees
  [trees specs]
  (map #(add-gamma-meta-to-tree % specs) trees))

(defn make-parental-tree
  "Given a species name which represents the 'hybrid' node, makes another tree
  that resolves the hybrid event."
  [trees spec]
  (reduce
   #(conj %1 %2 (permute-tree-at-spec %2 spec))
   '() 
   trees))

(defn make-parental-trees
  [tree specs]
  (let [p-trees (reduce #(make-parental-tree %1 %2) [tree] specs)]
    (add-gamma-meta-to-trees p-trees specs)))

(defn is-gamma-position-eq?
  "pos1 and pos2 are either 0,1,2, representing the position of the hybrid
  species.  Essentially, disregard (return true) if the spec is a hybrid (= 2)"
  [pos1 pos2]
  (or (= pos1 2) (= pos2 2) (= pos1 pos2)))

(defn is-parental-tree?
  "Given the gamma meta data for a parental tree, is the parental tree one that
  should be used to resolve hybridization in the hybrid tree given the hybrid tree's
  hybrid topology.  See gamma-topology doc-string for more info."
  [g-top g-meta]
  (let [logical-vals (map #(is-gamma-position-eq? (% g-top) (% g-meta)) (keys g-top))]
    (reduce #(and %1 %2) logical-vals)))

(defn p-trees-for-hybrid-tree
  "Uses the gamma meta data for each tree to find the appropriate parental
  trees to estimate the hybrid tree.  Uses the gamma topology map
  (e.g. {hyb1 0, hyb2 2}) to find these trees."
  [g-top p-trees]
  (filter #(is-parental-tree? (:gamma (meta %)) g-top) p-trees))

(defn gamma-top->num-hybrids
  [g-top]
  (reduce #(+ %1 (if-not (= %2 2) 1 0)) 0 (vals g-top)))

(defn gamma-topology->hybrid-tree-data
  "Given a gamma topology, i.e. the topology of a hybrid tree, find the gamma
  value that produces the mle."
  [gamma-top parental-trees gene-trees spec->lin theta]
  (let [n-hybrids (gamma-top->num-hybrids gamma-top)
        p-trees (p-trees-for-hybrid-tree gamma-top parental-trees)
        gamma-probs (find-gammas gene-trees p-trees spec->lin theta )
        k (compute-k n-hybrids (count (keys gamma-top)))
        aic (compute-aic (first gamma-probs) k)]
    (create-hybrid-tree-data (first p-trees) (n/vector-tree->newick-str (first p-trees))
                             (second gamma-probs) (first gamma-probs) k aic)))

(defn gamma-topology
  "For n hybridization events, there are 2^n parental trees, and
  3^n - 2^n hybrid trees.  Returns a seq of maps, where each map describes
  the topology of the hybrid tree.  For example, with n = 2 there are 4 parental
  trees and 5 hybrid trees.  Each hybrid tree has a unique topology described by
  a map of hybrid positions, e.g. {hyb1 0, hyb2 2}, which means that the first hybrid
  species is resolved to being a left child, and the 2 means that the second hybrid
  species needs to be estimated from the parental trees, namely that the two parental
  trees having a topology of {hyb1 0, hyb2 0} and {hyb1 0, hyb2 1} are used to
  estimate the gamma of the second hybrid event."
  [specs]
  (let [g-seq (filter #(>= (.indexOf % 2) 0) (c/selections [0,1,2] (count specs)))]
    (map #(apply hash-map (interleave specs %)) g-seq)))

(defn parental-trees->parental-data
  [parental-trees gene-trees spec->lin theta]
  (map (fn [p-tree]
            (let [lik (l/calc-lik gene-trees p-tree spec->lin theta)
                  newick (n/vector-tree->newick-str p-tree)
                  k (compute-k 0 (count (keys spec->lin)))
                  aic (compute-aic lik k)]
              (create-hybrid-tree-data p-tree newick nil nil lik k aic)
              ))
       parental-trees))
