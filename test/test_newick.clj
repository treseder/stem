(ns test-newick
  (:use [stem.newick])
  (:use [clojure.test]))


(def t-str "((one:1,two:2):0.2,five:5)")
(def bad-test-str "((one:1,two:2):0.2,)")
(def test-str "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7)" )
(def new-str "(((Name1:0.00123,Name2:0.00123):0.00123,(Name3:0.00121,Name4:0.00121) :0.00125):0.0010,((Name5:0.0010,Name6:0.0010):0.0014,(MyName7:0.0012,Name8:0.0012) :0.0012):0.00106)")


(deftest test-build-tree-from-newick-str
  (is (= (build-tree-from-newick-str t-str) [:root [:0.2 [two:2] [one:1]] [five:5]])))
