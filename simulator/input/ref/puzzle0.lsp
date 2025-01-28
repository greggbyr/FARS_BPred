;; $Header: puzzle.cl,v 1.1 87/12/30 09:46:06 layer Exp $
;; $Locker:  $

  (defconstant puzzle-size 511.)	
  (defconstant puzzle-classmax 3.)
  (defconstant puzzle-typemax 12.)

(defvar **iii** 0)
(defvar **kount** 0)
(defvar puzzle-d 8.)
(proclaim '(type fixnum **iii** **kount** puzzle-d))

(defvar piececount (make-array-init (1+ puzzle-classmax)  0))
(defvar puzzle-class (make-array-init (1+ puzzle-typemax) 0))
(defvar piecemax (make-array-init (1+ puzzle-typemax) 0))
(defvar puzzle (make-array (1+ puzzle-size)))
(defvar puzzle-p (make-array-2d (list (1+ puzzle-typemax) (1+ puzzle-size))))

(proclaim '(type simple-vector 
		 piececount puzzle-class piecemax puzzle))

		 
(defun fit (i j)
  (declare (type fixnum i j))
  (let ((end (aref piecemax i)))
    (declare (type fixnum end))
    (do ((k 0 (the fixnum (1+ k))))
	((> k end) t)
      (declare (type fixnum k))
      (cond ((aref (aref puzzle-p i) k)
	     (cond ((aref puzzle (the fixnum (+ j k)))
		    (return nil))))))))

(defun place (i j)
  (declare (type fixnum i j))
  (let ((end (aref piecemax i)))
    (declare (type fixnum end))
    (do ((k 0 (the fixnum (1+ k))))
	((> k end))
      (declare (type fixnum k))
      (cond ((aref (aref puzzle-p i) k) 
	     (setf (aref puzzle (the fixnum (+ j k))) t))))
    (setf (aref piececount (aref puzzle-class i)) 
	  (the fixnum
	       (- (the fixnum
		       (aref piececount (aref puzzle-class i))) 1)))
    (do ((k j (the fixnum (1+ k))))
	((> k puzzle-size)
	 (terpri)
	 (princ "Puzzle filled")
	 0)
      (declare (type fixnum k))
      (cond ((not (aref puzzle k))
	     (return k))))))


(defun puzzle-remove (i j)
  (declare (type fixnum i j))
  (let ((end (aref piecemax i)))
    (declare (type fixnum end))
    (do ((k 0 (the fixnum (1+ k))))
	((> k end))
      (declare (type fixnum k))
      (cond ((aref (aref puzzle-p i) k)
	     (setf (aref puzzle (the fixnum (+ j k)))  nil))))
    (setf (aref piececount (aref puzzle-class i))
	  (+ (the fixnum (aref piececount (aref puzzle-class i))) 1))))

(defun trial (j)
  (declare (type fixnum j))
  (let ((k 0))
    (declare (type fixnum k))
    (do ((i 0 (the fixnum (1+ i))))
	((> i puzzle-typemax) 
	 (setq **kount** (the fixnum (1+ **kount**))) nil)
      (declare (type fixnum i))
      (cond ((not (= (the fixnum (aref piececount (aref puzzle-class i))) 0))
	     (cond ((fit i j)
		    (setq k (place i j))
		    (cond ((or (trial k)
			       (= k 0))
			   (setq **kount** (the fixnum (+ **kount** 1)))
			   (return t))
			  (t (puzzle-remove i j))))))))))

(defun definepiece (iclass ii jj kk)
  (declare (type fixnum ii jj kk))
  (let ((index 0))
    (declare (type fixnum index))
    (do ((i 0 (the fixnum (1+ i))))
	((> i ii))
      (declare (type fixnum i))
      (do ((j 0 (the fixnum (1+ j))))
	  ((> j jj))
	(declare (type fixnum j))
	(do ((k 0 (the fixnum (1+ k))))
	    ((> k kk))
	  (declare (type fixnum k))
	  (setq index  
	    (+ i 
	       (the fixnum 
		    (* puzzle-d 
		       (the fixnum 
			    (+ j 
			       (the fixnum 
				    (* puzzle-d k))))))))
	  (setf (aref (aref puzzle-p **iii**) index)  t))))
    (setf (aref puzzle-class **iii**) iclass)
    (setf (aref piecemax **iii**) index) 
    (cond ((not (= **iii** puzzle-typemax))
	   (setq **iii** (the fixnum (+ **iii** 1)))))))

(defun puzzle-start ()
  (do ((m 0 (the fixnum (1+ m))))
      ((> m puzzle-size))
    (declare (type fixnum m))
    (setf (aref puzzle m) t))
  (do ((i 1 (the fixnum (1+ i))))
      ((> i 5))
    (declare (type fixnum i))
    (do ((j 1 (the fixnum (1+ j))))
	((> j 5))
      (declare (type fixnum j))
      (do ((k 1 (the fixnum (1+ k))))
	  ((> k 5))
	(declare (type fixnum k))
	(setf (aref puzzle 
		    (+ i 
		       (the fixnum
			    (* puzzle-d 
			       (the fixnum
				    (+ j 
				       (the fixnum
					    (* puzzle-d k))))))))
	      nil)))) 
  (do ((i 0 (the fixnum (1+ i))))
      ((> i puzzle-typemax))
    (declare (type fixnum i))
    (do ((m 0 (the fixnum (1+ m))))
	((> m puzzle-size))
      (declare (type fixnum m))
      (setf (aref (aref puzzle-p i) m)  nil)))
  (setq **iii** 0)
  (definepiece 0 3 1 0)
  (definepiece 0 1 0 3)
  (definepiece 0 0 3 1)
  (definepiece 0 1 3 0)
  (definepiece 0 3 0 1)
  (definepiece 0 0 1 3)
  
  (definepiece 1 2 0 0)
  (definepiece 1 0 2 0)
  (definepiece 1 0 0 2)
  
  (definepiece 2 1 1 0)
  (definepiece 2 1 0 1)
  (definepiece 2 0 1 1)
  
  (definepiece 3 1 1 1)
  
  (setf (aref piececount 0) 13.)
  (setf (aref piececount 1) 3)
  (setf (aref piececount 2) 1)
  (setf (aref piececount 3) 1)
  (setq **kount** 0)
  (let ((m (+ 1 (the fixnum (* puzzle-d (the fixnum (+ 1 puzzle-d))))))
	(n 0))
    ;;(declare (type fixnum m n **kount**))
    (cond ((fit 0 m) (setq n (place 0 m)))
	  (t (format t "~%Error.")))
    (cond ((trial n) 
	    (print (list 'Success 'in **kount**)))
	   ;;(format t "~%Success in:") (print **kount**))
	  (t (format t "~%Failure.")))))

(defun testpuzzle ()
  (time (puzzle-start)))

(testpuzzle)
(print '(done puzzle))
