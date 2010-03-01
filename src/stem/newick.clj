(ns stem.newick
  (:require [clojure.contrib.str-utils2 :as s]
            [stem.util :as util]))

(defn build-tree [el]
  (if (= (count el) 1)
    [(first el)]
    [(nth el 2)
     (build-tree (nth el 1))
     (build-tree (nth el 0))]))

(defn prep-newick-str
  "To make the tree easier for parsing the commas are replaced with ')('."
  [n-str]
  (str "(" (s/replace n-str "," ")(") ")"))

(defn build-tree-from-newick-str
  "Parses s and builds a binary tree structure as a vector of
  vectors.  s must conform to proper newick format.  If the
  root is named, it is ignored.  Leafs are assumed to contain a
  name:number, where number is optional.  All interior nodes must
  have a number, but need not be named."
  [s]
  (let [prepped-str (prep-newick-str s)
        lst (read-string prepped-str)]
    (try
     [:root (build-tree (first lst)) (build-tree (second lst))]
     (catch Exception _
       (util/abort "An error occured parsing the newick tree.")))))
