(ns stem.gene-tree
  (:require [stem.util :as util]
            [stem.newick :as newick]
            [clojure.contrib.str-utils2 :as s]))

(defstruct gene-tree :vec-tree :rate)

(defn create-gene-tree [vec-tree rate]
  (struct gene-tree vec-tree rate))

(defn parse-rate-and-tree [s]
  (re-find #"\[(\d*\.*\d)\](.*)" s))


(defn create-struct-from-parts
  [[_ rate-str newick-str]]
  (create-gene-tree
   (newick/build-tree-from-newick-str newick-str)
   (Double/parseDouble rate-str)))

(defn parse-gene-tree-str [s]
  (->> s
       (parse-rate-and-tree)
       (create-struct-from-parts)))

(defn parse-gene-tree-file [file-name]
  (let [file-str (util/remove-whitespace (slurp file-name))
        coll-newick-strs (s/split file-str #";")]
    (map parse-gene-tree-str coll-newick-strs)))
