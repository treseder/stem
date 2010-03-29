(ns stem.gene-tree
  (:require [stem.util :as util]
            [stem.newick :as newick]
            [clojure.contrib.str-utils2 :as s]))

(defstruct gene-tree :vec-tree :rate)

(defn create-gene-tree [vec-tree rate]
  (struct gene-tree vec-tree rate))

(defn parse-rate-and-tree [s]
  (re-find #"(\[(\d*\.*\d)\])?(.*)" s))


(defn create-struct-from-parts
  [rate-str newick-str]
  (create-gene-tree
   (newick/build-tree-from-newick-str newick-str)
   (Double/parseDouble rate-str)))

(defn parse-gene-tree-str
  ([s & r]
  (let [[_ _ rate n-str] (parse-rate-and-tree)]
    (create-struct-from-parts (or r rate) n-str))))

(defn parse-gene-tree-file
  "Returns a sequence of gene-tree structs found in file-name"
  [file-name & rate]
  (let [file-str (util/remove-whitespace (slurp file-name))
        newick-strs (s/split file-str #";")]
    (map #(parse-gene-tree-str % rate) newick-strs)))

(defn get-gene-trees [file-map]
  (if file-map
    (reduce (fn [s [file-name rate]] (concat s (parse-gene-tree-file file-name rate))) [] file-map)
    ;assumes all genetrees are in one file called genetrees.tre
    (parse-gene-tree-file "genetrees.tre")))
