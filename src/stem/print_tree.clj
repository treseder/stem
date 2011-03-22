(ns stem.print-tree)

;; Taken from Conrad Barski's Vijual library

(defn out 
  "like 'print' but without inserted spaces and fully flushed"
  [& more]
  (print (apply str more))
  (flush)
  (first more))

(defn fill 
  "Returns a string with the character repeated"
  ([c n]
     (apply str (take n (repeat c))))
  ([c]
     (fill c 1)))

(defn sp 
  "Returns 1->n spaces in a string"
  ([n]
     (out (fill \space n)))
  ([]
     (sp 1)))

(defn label-text 
  "Returns a text string representing the label for the item. It handles the special case of keywords, so that :oak-tree ==> 'oak tree'"
  [text]
     (if (keyword? text)
       (apply str (replace {\- \space} (name text)))
       (str text)))

(defn line-info-btree 
  "Takes 2 rows of nodes in the binary tree and figures out the arrangement of lines between the two rows."
  [top bottom]
  (letfn [[f [top tx bottom bx in-line x]
           (if (and (seq top) (or (empty? bottom) (< (+ tx (first (first top))) (+ bx (first (first bottom))))))
             (let [[ind wid val left right] (first top)]
                  (concat (if left
                            [[:lbottom] [:line (- (+ tx ind) bx 2)] [:ltop] [:space wid] [:nop]] 
                            [[:space (- (+ tx ind wid) x)] [:nop]])
                          (f (rest top) (+ tx ind wid) bottom bx right (+ tx ind wid))))
             (when (seq bottom)
               (let [[ind wid val left right] (first bottom)]
                   (concat (if in-line
                             [[:rtop] [:line (- (+ bx ind) tx 2)] [:rbottom] [:space wid] [:nop]]
                             [[:space (- (+ bx ind wid) x)] [:nop]])
                           (f top tx (rest bottom) (+ bx ind wid) false (+ bx ind wid))))))]]
    (f top 0 bottom 0 false 0)))

(defn btree-row-wid [row]
  "Figures out the width of a row of the a binary tree."
  (apply + (map (fn [[a b c]]
                  (+ a b))
                row)))

(defn layout-btree 
  "Takes a binary tree and converts it into rows, with nodes at the same level of the tree in the same row. The data for each item in each row is a vector, formatted as [indentation width text left-children? right-children?"
  [btree]
  (if btree
    (let [[cur no yes] btree
          lno (layout-btree no)
          lyes (layout-btree yes)
          wno (apply max 0 (map btree-row-wid lno))
          wyes (apply max 0 (map btree-row-wid lyes))
          wid (+ (count (label-text cur)) 4)
          node-off (if (empty? lno)
                     0
                     (+ (btree-row-wid (first lno)) 2))
          yes-off (max (inc wno) (+ node-off wid 2))]
      (cons [[node-off wid cur (not (empty? lno)) (not (empty? lyes))]]
            (let [m (max (count lno) (count lyes))]
              (map (fn [rno ryes]
                     (if ryes
                       (let [[[a b c d e] & t] ryes]
                           (concat rno (cons [(- (+ a yes-off) (btree-row-wid rno)) b c d e] t)))
                       rno))
                   (take m (concat lno (repeat nil)))
                   (take m (concat lyes (repeat nil)))))))
    []))

(defn render-btree [rows]
  "Renders a binary tree. The rows handed to it contain the indentation already, so this function is concerned solely with rendering this row data to ASCII"
  (let [x rows]
    (loop [rows x]
      (when (seq rows)
        (let [row (first rows)]
          (doseq [[ind w str] row]
            (out (fill \space ind) \+ (fill \- (- w 2)) \+))
          (newline)
          (doseq [[ind w str] row]
            (dotimes [i ind]
              (out \space))
            (out "| " (label-text str) " |"))
          (newline)
          (doseq [[ind w str] row]
            (out (fill \space ind) \+ (fill \- (- w 2)) \+))
          (newline)
          (when (seq (rest rows))
            (let [li (line-info-btree row (first (rest rows)))]
              (doseq [[type n] li]
                (condp = type
                      :space (sp n)
                      :lbottom (sp)
                      :line (out (fill \_ n))
                      :ltop (out (fill \/))
                      :rtop (out (fill \\))
                      :rbottom (sp)
                      :nop nil))
              (newline)
              (doseq [[type n] li]
                (condp = type
                      :space (sp n)
                      :lbottom (out (fill \/))
                      :line (sp n)
                      :ltop (sp)
                      :rtop (sp)
                      :rbottom (out (fill \\))
                      :nop nil))
              (newline))))
        (recur (rest rows))))))

(defn draw-binary-tree
  "Draws a binary tree to the console. Nodes are in the form [text left right] where 'left' and 'right' are optional fields containing the children of this node."
  [btree]
  (render-btree (layout-btree btree)))
