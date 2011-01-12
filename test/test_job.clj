(ns test-job
  (:use [stem.job] [clojure.test]))

(defn get-results [settings]
  (:results (-> (create-job settings) (pre-run-check) (run))))

(deftest test-simple-mle
  (testing "simple mle and ml tree"
    (let [results (get-results "simplesettings.yaml")]
      (is (= -52.43701947216076 (:mle results)))
      (is (= "(Species1:1.20000,(Species4:1.10000,(Species2:1.00000,Species3:1.00000):0.10000):0.10000);" (:species-tree results))))))


(comment
  (deftest test-tied-trees
   (testing "tied trees"
     (let [results (get-results "settings-ties.yaml")]
       (is (= 3 (count (:tied-trees results))))))))

(deftest test-missing-data
  (testing "missing data"
    (let [results (get-results "settings-dan.yaml")]
      (is (= -16089.0014171403826 (:mle results)))
      (is (= "(midgroup2:31.27000,(outgroup:8.06000,(midgroup:7.57000,ingroup:7.57000):0.49000):23.21000);" (:species-tree results))))))

(deftest test-search
  (testing "search with 100 iterations, seed = 1"
    (let [results (get-results "settings-rattle.yaml")]
      (is (= -10786.177600241292 ((first (:best-trees results)) 0)))
      (is (= -13609.008520241292 ((last (:best-trees results)) 0)))      
      (is (= 10 (count (:best-trees results)))))))

(deftest user-trees
  (testing "user trees"
    (let [results (get-results "simplesettings-user.yaml")]
      (is (= '(-1610.6688720371628 -1597.5895120371629 -56.63839947216076)
             (:user-liks results)))
      (is (= '(-53.63717947216077 -56.63683947216076 -56.63719947216075)
             (:optim-liks results))))))

(def user-test-trees
     "(Species1:1.40000,(Species3:1.10000,(Species2:0.5,Species4:0.5):0.6):0.30000);
((Species2:1.00000,Species3:1.00000):0.50000,(Species4:1.09995,Species1:1.09995):0.50005);
(((Species2:1.00000,Species3:1.00000):0.09995,Species1:1.09995):0.00005,Species4:1.10000);")
