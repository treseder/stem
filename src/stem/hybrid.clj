(ns stem.hybrid
  (:use [stem.constants] [clojure.set])
  (:require [stem.newick :as n]
            [stem.lik :as l]
            [stem.util :as u]
            [clojure.zip :as z]))


(defn compute-aic
  [lik k]
  (+ (* (- 2) lik) (* 2 k)))

(defn compute-k
  [hybs specs]
  (+ hybs (- specs 1)))

(defn find-gamma-2
  [g-trees spec-trees spec-to-lin theta])

(defn reduce-gene-trees
  [lik-memo g-trees spec-trees spec-to-lin theta gamma]
  (let [[spt1 spt2] spec-trees
        prob (reduce
              (fn [total gene-tree]
                (let [p1 (* gamma (Math/exp (lik-memo (:vec-tree gene-tree) spt2 spec-to-lin (/ 2 theta))))
                      p2 (* (- 1 gamma) (Math/exp (lik-memo (:vec-tree gene-tree) spt1 spec-to-lin (/ 2 theta))))]
                  ;(println (str "p1 p2: " p1 " " p2))
                  (+ total (Math/log (+ p1 p2)))))
              0 g-trees)]
    prob))

(defn find-gamma-1
  [g-trees spec-trees spec-to-lin theta]
  (let [lik-memo (memoize l/calc-lik-for-tree)
        gamma-fn (partial reduce-gene-trees lik-memo g-trees spec-trees spec-to-lin theta)]
    (reduce
     (fn [[curr-lik gamma] new-gamma]
       (let [new-lik (gamma-fn new-gamma)]
         (if (> new-lik curr-lik)
           [new-lik [new-gamma]]
           [curr-lik gamma])))
     [Double/NEGATIVE_INFINITY [0]]
     (range 0 1 0.01))))


(defn find-gammas
  [g-trees spec-trees spec-to-lin theta]
  (if (= (count spec-trees) 2)
    (find-gamma-1 g-trees spec-trees spec-to-lin theta)
    (find-gamma-2 g-trees spec-trees spec-to-lin theta)))


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
