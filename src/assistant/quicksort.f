      recursive subroutine quicksort(a, v, first, last)
        implicit none
        integer  a(*), i1, v(*), x, t
        integer first, last
        integer i, j
      
        x = v( (first+last) / 2 )
        i = first
        j = last
        do
           do while (v(i) < x)
              i=i+1
           end do
           do while (x < v(j))
              j=j-1
           end do
           if (i >= j) exit
           t  = v(i);  v(i) = v(j);  v(j) = t; 
           i1 = a(i) ; a(i) = a(j) ; a(j) = i1;
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, v, first, i-1)
        if (j+1 < last)  call quicksort(a, v, j+1, last)
      end subroutine quicksort
