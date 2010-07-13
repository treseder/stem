(ns test-newick
  (:use [stem.newick])
  (:use [clojure.test]))


(def newick-1 "((one:1,two:2):0.2,five:5)")
(def tree-1 [{:name "stem-root", :time 0.0, :c-time 5.0, :desc #{"one" "two" "five"}} [{:name "internal", :time 0.2, :c-time 2.0, :desc #{"one" "two"}} [{:name "one", :time 1, :c-time 0.0, :desc #{}}] [{:name "two", :time 2, :c-time 0.0, :desc #{}}]] [{:name "five", :time 5, :c-time 0.0, :desc #{}}]])

(def newick-2 "(((Name1:0.00123,Name2:0.00123):0.00123,(Name3:0.00121,Name4:0.00121):0.00125):0.0010,((Name5:0.0010,Name6:0.0010):0.0014,(MyName7:0.0012,Name8:0.0012):0.0012):0.00106)")
(def tree-2 [{:name "stem-root", :time 0.0, :c-time 0.0034600000000000004, :desc #{"MyName7" "Name1" "Name2" "Name3" "Name4" "Name5" "Name6" "Name8"}} [{:name "internal", :time 0.0010, :c-time 0.00246, :desc #{"Name1" "Name2" "Name3" "Name4"}} [{:name "internal", :time 0.00123, :c-time 0.00123, :desc #{"Name1" "Name2"}} [{:name "Name1", :time 0.00123, :c-time 0.0, :desc #{}}] [{:name "Name2", :time 0.00123, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.00125, :c-time 0.00121, :desc #{"Name3" "Name4"}} [{:name "Name3", :time 0.00121, :c-time 0.0, :desc #{}}] [{:name "Name4", :time 0.00121, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.00106, :c-time 0.0024000000000000002, :desc #{"MyName7" "Name5" "Name6" "Name8"}} [{:name "internal", :time 0.0014, :c-time 0.0010, :desc #{"Name5" "Name6"}} [{:name "Name5", :time 0.0010, :c-time 0.0, :desc #{}}] [{:name "Name6", :time 0.0010, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.0012, :c-time 0.0012, :desc #{"MyName7" "Name8"}} [{:name "MyName7", :time 0.0012, :c-time 0.0, :desc #{}}] [{:name "Name8", :time 0.0012, :c-time 0.0, :desc #{}}]]]])

(def newick-3 "((((Name1:0.00123,Name2:0.00123):0.00133,(Name3:0.0012,Name4:0.0012):0.00134):0.0003,(Name5:0.0010,Name6:0.0010):0.00186):0.00064,(MyName7:0.0011,Name8:0.0011):0.0024)")
(def tree-3 [{:name "stem-root", :time 0.0, :c-time 0.0035, :desc #{"MyName7" "Name1" "Name2" "Name3" "Name4" "Name5" "Name6" "Name8"}} [{:name "internal", :time 6.4E-4, :c-time 0.00286, :desc #{"Name1" "Name2" "Name3" "Name4" "Name5" "Name6"}} [{:name "internal", :time 3.0E-4, :c-time 0.0025599999999999998, :desc #{"Name1" "Name2" "Name3" "Name4"}} [{:name "internal", :time 0.00133, :c-time 0.00123, :desc #{"Name1" "Name2"}} [{:name "Name1", :time 0.00123, :c-time 0.0, :desc #{}}] [{:name "Name2", :time 0.00123, :c-time 0.0, :desc #{}}]] [{:name "internal", :time 0.00134, :c-time 0.0012, :desc #{"Name3" "Name4"}} [{:name "Name3", :time 0.0012, :c-time 0.0, :desc #{}}] [{:name "Name4", :time 0.0012, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.00186, :c-time 0.0010, :desc #{"Name5" "Name6"}} [{:name "Name5", :time 0.0010, :c-time 0.0, :desc #{}}] [{:name "Name6", :time 0.0010, :c-time 0.0, :desc #{}}]]] [{:name "internal", :time 0.0024, :c-time 0.0011, :desc #{"MyName7" "Name8"}} [{:name "MyName7", :time 0.0011, :c-time 0.0, :desc #{}}] [{:name "Name8", :time 0.0011, :c-time 0.0, :desc #{}}]]])

(def ntree-6 "(Species6:6.303000,((Species1:3.052820,Species5:3.052820):1.143503,(Species2:2.023739,(Species3:1.117750,Species4:1.117750):0.905989):2.172584):2.106677)")
(def ntree-6-mle -56.463282)

(def ntree-12 "((Species6:6.068630,((Species1:3.052820,Species5:3.052820):0.983191,(Species2:2.023739,(Species3:1.114170,Species4:1.114170):0.909569):2.012272):2.032619):2.041061,(Species7:6.032850,((Species8:3.502780,Species9:3.502780):0.512300,(Species10:2.051490,(Species11:1.304710,Species12:1.304710):0.746780):1.963590):2.017770):2.076841)")
(def tree-12-mle -115.524690)

(def ntree-18 "(Species6:6.068630,(Species7:3.100870,(((Species1:2.018190,Species2:2.018190):0.005549,(Species3:1.015080,Species4:1.015080):1.008659):0.087421,(Species5:1.090940,Species8:1.090940):1.020220):0.989710):2.967760)")
(def ntree-18-mle -4211.707401)

(deftest test-build-tree-from-newick-str
  (is (= (build-tree-from-newick-str newick-1 1 1) tree-1))
  (is (= (build-tree-from-newick-str newick-2 1 1) tree-2))
  (is (= (build-tree-from-newick-str newick-3 1 1) tree-3)))


