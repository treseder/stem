(ns test-newick
  (:use [stem.newick])
  (:use [clojure.test]))


(def newick-1 "((one:1,two:2):0.2,five:5)")
(def tree-1 [{:name "stem-root", :time 0.0, :c-time 5.0, :desc #{"one" "two" "five"}} [{:name "internal", :time 0.2, :c-time 2.0, :desc #{"one" "two"}} [{:name "one", :time 1, :c-time 0.0, :desc #{}}] [{:name "two", :time 2, :c-time 0.0, :desc #{}}]] [{:name "five", :time 5, :c-time 0.0, :desc #{}}]])

(def newick-2 "(((Name1:0.00123,Name2:0.00123):0.00123,(Name3:0.00121,Name4:0.00121):0.00125):0.0010,((Name5:0.0010,Name6:0.0010):0.0014,(MyName7:0.0012,Name8:0.0012):0.0012):0.00106)")
(def tree-2 [{:name "stem-root", :time 0.0, :c-time 0.0034600000000000004, :desc #{"MyName7" "Name1" "Name2" "Name3" "Name4" "Name5" "Name6" "Name8"}} [{:name "internal", :time 0.0010, :c-time 0.00246, :desc #{"Name1" "Name2" "Name3" "Name4"}} [{:name "internal", :time 0.00123, :c-time 0.00123, :desc #{"Name1" "Name2"}} [{:name "Name1", :time 0.00123, :c-time 0.0, :desc #{}}] [{:name "Name2", :time 0.00123, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.00125, :c-time 0.00121, :desc #{"Name3" "Name4"}} [{:name "Name3", :time 0.00121, :c-time 0.0, :desc #{}}] [{:name "Name4", :time 0.00121, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.00106, :c-time 0.0024000000000000002, :desc #{"MyName7" "Name5" "Name6" "Name8"}} [{:name "internal", :time 0.0014, :c-time 0.0010, :desc #{"Name5" "Name6"}} [{:name "Name5", :time 0.0010, :c-time 0.0, :desc #{}}] [{:name "Name6", :time 0.0010, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.0012, :c-time 0.0012, :desc #{"MyName7" "Name8"}} [{:name "MyName7", :time 0.0012, :c-time 0.0, :desc #{}}] [{:name "Name8", :time 0.0012, :c-time 0.0, :desc #{}}]]]])

(def newick-3 "((((Name1:0.00123,Name2:0.00123):0.00133,(Name3:0.0012,Name4:0.0012):0.00134):0.0003,(Name5:0.0010,Name6:0.0010):0.00186):0.00064,(MyName7:0.0011,Name8:0.0011):0.0024)")
(def tree-3 [{:name "stem-root", :time 0.0, :c-time 0.0035, :desc #{"MyName7" "Name1" "Name2" "Name3" "Name4" "Name5" "Name6" "Name8"}} [{:name "internal", :time 6.4E-4, :c-time 0.00286, :desc #{"Name1" "Name2" "Name3" "Name4" "Name5" "Name6"}} [{:name "internal", :time 3.0E-4, :c-time 0.0025599999999999998, :desc #{"Name1" "Name2" "Name3" "Name4"}} [{:name "internal", :time 0.00133, :c-time 0.00123, :desc #{"Name1" "Name2"}} [{:name "Name1", :time 0.00123, :c-time 0.0, :desc #{}}] [{:name "Name2", :time 0.00123, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.00134, :c-time 0.0012, :desc #{"Name3" "Name4"}} [{:name "Name3", :time 0.0012, :c-time 0.0, :desc #{}}] [{:name "Name4", :time 0.0012, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.00186, :c-time 0.0010, :desc #{"Name5" "Name6"}} [{:name "Name5", :time 0.0010, :c-time 0.0, :desc #{}}] [{:name "Name6", :time 0.0010, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.0024, :c-time 0.0011, :desc #{"MyName7" "Name8"}} [{:name "MyName7", :time 0.0011, :c-time 0.0, :desc #{}}] [{:name "Name8", :time 0.0011, :c-time 0.0, :desc #{}}]]])


(deftest test-build-tree-from-newick-str
  (is (= (build-tree-from-newick-str newick-1 1 1) tree-1))
  (is (= (build-tree-from-newick-str newick-2 1 1) tree-2))
  (is (= (build-tree-from-newick-str newick-3 1 1) tree-3)))
