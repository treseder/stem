(ns stem.lik-tree
  (:require [clojure.set :as set]
            [stem.util :as u]
            [stem.newick :as newick]
            [stem.messages :as m]
            [clojure.contrib.combinatorics :as comb]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;functions for building the coalescent matrix ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn get-sorted-name-by-index [idx-map]
  (sort-by #(idx-map %) (keys idx-map)))

(defn name->indexes
  "Takes a seq of lineage names and turns them into the indexes of the matrix"
  [lins index-map]
  (map index-map lins))

(defn fill-matrix-for-lineage
  [i-idx j-idxs matrix c-time]
  (doseq [j-idx j-idxs]
    (if (> i-idx j-idx) ;; matrix is symmetric - only fill half
      (u/aset! matrix j-idx i-idx c-time)
      (u/aset! matrix i-idx j-idx c-time))))

(defn check-lin [idx name]
  (if (nil? idx) (u/abort (str (m/e-strs :missing-lin) name) nil))
  idx)

(defn rec-fill-time-matrix
  [node matrix lin-index parent-desc-set parent-c-time]
  (let [[{name :name  desc-set :desc c-time :c-time} l-node r-node] node
        ;; if leaf consider self the only desc
        checked-desc-set (if (empty? desc-set) #{name} desc-set)
        diff-set (set/difference parent-desc-set checked-desc-set)
        js (name->indexes diff-set lin-index)]
    (doseq [i checked-desc-set]
      (fill-matrix-for-lineage (check-lin (lin-index i) i) js matrix parent-c-time))
    (if-not (newick/is-leaf? node)
      (do (rec-fill-time-matrix l-node matrix lin-index desc-set c-time)
          (rec-fill-time-matrix r-node matrix lin-index desc-set c-time)))))

(defn fill-time-matrix-for-tree
  "Returns a matrix filled with coalescent times for each of the lineages.
  The matrix is upper triangular."
  [tree matrix lin-index]
  (rec-fill-time-matrix tree matrix lin-index (:desc (first tree)) 0)
  matrix)

(defn gene-trees-to-matrices
  "Takes a seq of gene-trees and returns a seq of matrices"
  [trees lin-indexes]
  (let [len (count lin-indexes)]
    (map #(fill-time-matrix-for-tree (:vec-tree %)
                                     (u/make-stem-array len len)
                                     lin-indexes)
         trees)))

(defn reduce-matrices
  "Assumes matrices are pxp, upper-triangular.  Checks each cell in each
  matrix and stores the lesser in a.  Returns a."
  [a b]
  (let [p (count a)]
    (dorun
     (for [i (range p) j (range (+ i 1) p)]
       (let [a-val (u/aget! a i j)
             b-val (u/aget! b i j)]
         (if (or (zero? a-val) (and (< b-val a-val) (u/not-zero? b-val)))
           (u/aset! a i j b-val)))))
    a))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions for generating species matrix;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn fill-matrix-for-species [i j c-time s-matrix]
  (let [current-val (if (> i j) (u/aget! s-matrix j i)  (u/aget! s-matrix i j) )]
    (if (and (u/not-zero? c-time) (or (< c-time current-val) (zero? current-val)))
      (if (< i j)
        (u/aset! s-matrix i j c-time)
        (u/aset! s-matrix j i c-time)))))

(defn to-spec-matrix
  [l-matrix {:keys [lin-to-index spec-set spec-to-index lin-to-spec]}]
  (let [s-count (count spec-set)
        s-matrix (u/make-stem-array s-count s-count)
        l-count (count lin-to-spec)
        index-to-lin (zipmap (vals lin-to-index) (keys lin-to-index))]
    ;; for loop only traverses the upper right triangle of the matrix
    (dorun (for [i (range l-count) j (range (+ i 1) l-count)]
       (let [i-lin (index-to-lin i)
             j-lin (index-to-lin j)
             i-spec (spec-to-index (lin-to-spec i-lin))
             j-spec (spec-to-index (lin-to-spec j-lin))]
         (if-not (= i-spec j-spec)
           (fill-matrix-for-species i-spec j-spec (u/aget! l-matrix i j) s-matrix)))))
    s-matrix))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions for generating newick str from species matrix;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn sort-vec-of-vec [coll]
  (sort #(compare (first %1) (first %2)) coll))


(defn matrix->sorted-list
  "Returns a list sorted by elements of the matrix.
  Assumes matrix is upper triangular where diag elements are zero"
  [m idx-to-name]
  (let [size (count m)
        v (transient [])]
    (dorun
     (for [i (range size) j (range (+ i 1) size)]
       (let [time (u/aget! m i j)]
         (conj! v [time (idx-to-name i) (idx-to-name j)]))))
    (sort #(compare (first %1) (first %2)) (persistent! v))))

(defn partial-find
  "Each key in m is a set of keys. Return value of map entry where k
  is a member of such set."
  [k m]
  (first (filter (fn [[ks v]] (contains? ks k)) m)))

(defn nil->node [n v]
  (if (vector? v) v [n]))

(defn nil->set [n aset]
  (if (nil? aset) #{n} aset))

(defn find-and-merge-nodes
  "Given a node and a map of nodes merge the node with its appropriate
  ancestor.
  node-map is of the form {#{s1 s2} [.55 [s1][s2]]}
  node is a vector of size 3 that contains the coalescent time of the two species:
  [0.55 s1 s2].  Ancestors are all internal nodes, with the species being leaves
  of the tree."
  [node-map node]
  (let [[time l-name r-name] node
        [l-descs l-tree] (partial-find l-name node-map)
        [r-descs r-tree] (partial-find r-name node-map)]
    (if-not (or (contains? l-descs r-name) (contains? r-descs l-name))
      (let [comb-tree (if (nil? l-tree)
                        [time (nil->node l-name l-tree) (nil->node r-name r-tree)]
                        [time (nil->node r-name r-tree) (nil->node l-name l-tree)])
            new-descs (set/union (nil->set l-name l-descs) (nil->set r-name r-descs))]
        (-> node-map (dissoc l-descs r-descs) (assoc new-descs comb-tree)))
      node-map)))

(defn tree-from-seq [s]
  (reduce #(find-and-merge-nodes %1 %2) {} s))

(defn remove-unneeded-nodes
  [spec-set spec-lst]
  (filter (fn [[time l r]]
            (not (and (contains? spec-set l)
                      (contains? spec-set r))))
          spec-lst))

(defn get-specs-from-seq
  "Sequence is of the form:
  ([time spec-name1 spec-name2] [etc...])"
  [lst]
  (reduce (fn [spec-set [_ l-name r-name]]
            (conj spec-set l-name r-name))
          #{}
          lst))

(defn remove-unneeded-species-nodes
  [spec-lst]
  (reduce (fn [[spec-set cum-lst] lst]
            (let [filt-lst (remove-unneeded-nodes spec-set lst)
                  filt-set (get-specs-from-seq filt-lst)]
              (if-not (empty? filt-lst)
                [(set/union spec-set filt-set) (conj cum-lst filt-lst)]
                ; keeps at least one from each group
                [spec-set (conj cum-lst (list (first lst)))])))
          [#{} '()]
          spec-lst))

(defn get-list-permutations
  "The list of nodes used to build the ml species tree contains ties in how
  the tree should be built, with some nodes having the same coalescent time.
  Find all permutations of this list, thus finding all tied likelihood trees"
  [lst]
  (let [[_ clean-lst] (remove-unneeded-species-nodes
                       (vals (into (sorted-map) (group-by first lst))))]
    (map
     #(apply concat %)
     (apply comb/cartesian-product (map comb/permutations clean-lst)))))


(defn find-tied-trees [spec-lst]
  (->> spec-lst
       (get-list-permutations)
       (map #(second (first (tree-from-seq (sort-by first %1)))))
       (set)))

(defn get-lik-tree-parts
  ([gene-trees env] (get-lik-tree-parts gene-trees env true))
  ([gene-trees env find-tied-trees?]
   (let [gene-matrices  (gene-trees-to-matrices gene-trees (env :lin-to-index))
         min-gene-matrix (reduce reduce-matrices
                                 (u/make-stem-array (env :mat-size)) gene-matrices)
         spec-matrix (to-spec-matrix min-gene-matrix env)
         spec-lst (matrix->sorted-list spec-matrix (env :index-to-spec))
         tied-trees (if find-tied-trees? (find-tied-trees spec-lst) spec-lst)]
     {:min-gene-matrix min-gene-matrix
      :spec-matrix spec-matrix
      :tied-trees tied-trees})))

