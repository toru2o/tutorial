/*

start codon = ATG

*/
package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
	"math"
	"math/big"
	"net/http"
	"os"
	"strconv"
	"strings"
	"time"
)

// node番号は0, 1, 2, ...   1始まりにするとプログラムミスをしやすいので
func isConnected(adjacent [][]int, start, goal int) bool {
	path := []int{} //path := make([]int, 0, Size)
	ret := false

	member := func(n int, xs []int) bool {
		for _, x := range xs {
			if n == x {
				return true
			}
		}
		return false
	}

	var dfs func(int, []int)
	dfs = func(goal int, path []int) {
		if !ret {
			n := path[len(path)-1]
			if n == goal {
				ret = true
			} else {
				// node番号を1始まりとする場合は range adjacent[n-1]とする
				for _, x := range adjacent[n] {
					if !member(x, path) {
						dfs(goal, append(path, x))
					}
				}
			}
		}
	}

	dfs(goal, append(path, start))
	return ret

}

func convStoI(s string) int {
	n, err := strconv.Atoi(s)
	if err != nil {
		panic("err")
	}
	return n
}

func freadLines(fname string) []string {
	file, err := os.Open(fname)
	if err != nil {
		//os.Exit(1)
		panic("該当するファイルが見つかりません!")
	}
	defer file.Close()

	lines := []string{}
	sc := bufio.NewScanner(file)
	for i := 1; sc.Scan(); i++ {
		if err := sc.Err(); err != nil {
			fmt.Println(err)
			break
		}
		//fmt.Printf("%4d行目: %s\n", i, sc.Text())
		lines = append(lines, sc.Text())
	}
	return lines
}

func getIdData(fname string) []string {
	lt := freadLines(fname)
	lines := []string{}
	ss := ""
	for _, v := range lt {
		if v[0] == '>' {
			if ss != "" {
				lines = append(lines, ss)
				ss = ""
			}
			lines = append(lines, v)
		} else {
			ss += v
		}
	}
	lines = append(lines, ss)
	return lines
}

func fwriteLine(fname string, line string) {
	fout, err := os.Create(fname)
	if err != nil {
		panic("ファイル書き込みエラー")
	}
	defer fout.Close()

	fout.Write(([]byte)(line + "\r\n"))
}

func fwriteLines(fname string, lines []string) {
	fout, err := os.Create(fname)
	if err != nil {
		//os.Exit(1)
		panic("ファイル書き込みエラー")
	}
	defer fout.Close()

	for i := 0; i < len(lines); i++ {
		fout.Write(([]byte)(lines[i] + "\r\n"))
	}
	//fout.Close()
}

func MapString(f func(string) string, vs []string) []string {
	vsm := make([]string, len(vs))
	for i, v := range vs {
		vsm[i] = f(v)
	}
	return vsm
}

func combi(ss []int, m int) [][]int {
	seq := [][]int{}

	var combiSub func(xs []int)
	combiSub = func(xs []int) {
		if len(xs) == m {
			lt := []int{}
			for _, v := range xs {
				lt = append(lt, ss[v])
			}
			seq = append(seq, lt)
		} else {
			si := 0
			if len(xs) > 0 {
				si = xs[len(xs)-1] + 1
			}
			for i := si; i < len(ss); i++ {
				if len(xs)+len(ss)-i < m {
					break
				}
				combiSub(append(xs, i))
			}
		}
	}

	combiSub(make([]int, 0, m))
	return seq
}

func combinate(ss []string, m int) [][]string {
	seq := [][]string{}

	var combiSub func(xs []int)
	combiSub = func(xs []int) {
		if len(xs) == m {
			lt := []string{}
			for _, v := range xs {
				lt = append(lt, ss[v])
			}
			seq = append(seq, lt)
		} else {
			si := 0
			if len(xs) > 0 {
				si = xs[len(xs)-1] + 1
			}
			for i := si; i < len(ss); i++ {
				if len(xs)+len(ss)-i < m {
					break
				}
				combiSub(append(xs, i))
			}
		}
	}
	combiSub(make([]int, 0, m))
	return seq
}

// 文字列（スライスに格納）から重複を許してm個を取り出して並べる順列を出力
func permuteDup(ss []string, m int) [][]string {
	seq := [][]string{}

	var permSub func(xs []string)
	permSub = func(xs []string) {
		if len(xs) == m {
			lt := []string{}
			for _, v := range xs {
				lt = append(lt, v)
			}
			seq = append(seq, lt)
		} else {
			for i := 0; i < len(ss); i++ {
				permSub(append(xs, ss[i]))
			}
		}
	}

	permSub(make([]string, 0, m))
	return seq
}

// 文字列（スライスに格納）からm個を取り出して並べる順列を出力
func permutateStr(ss []string, m int) [][]string {
	seq := [][]string{}
	//used := make([]int, len(ss))
	member := func(a string, xs []string) bool {
		for i := 0; i < len(xs); i++ {
			if xs[i] == a {
				return true
			}
		}
		return false
	}

	var permSub func(xs []string)
	permSub = func(xs []string) {
		if len(xs) == m {
			lt := []string{}
			for _, v := range xs {
				lt = append(lt, v)
			}
			seq = append(seq, lt)
		} else {
			for i := 0; i < len(ss); i++ {
				if !member(ss[i], xs) {
					permSub(append(xs, ss[i]))
				}
			}
		}
	}

	permSub(make([]string, 0, m))
	return seq
}

// 整数列の順列をスライスで出力
func permutation(n, m int) [][]int {
	seq := [][]int{}

	member := func(n int, xs []int) bool {
		for _, x := range xs {
			if n == x {
				return true
			}
		}
		return false
	}

	var permSub func(n, m int, xs []int)
	permSub = func(n, m int, xs []int) {
		if len(xs) == m {
			//f(xs)
			lt := []int{}
			for _, v := range xs {
				lt = append(lt, v)
			}
			seq = append(seq, lt)
		} else {
			for i := 1; i <= n; i++ {
				if !member(i, xs) {
					permSub(n, m, append(xs, i))
				}
			}
		}
	}

	permSub(n, m, make([]int, 0, m))
	return seq
}

func getBody(url string) string {
	//url := "http://www.uniprot.org/uniprot/B5ZC00.fasta"
	res, err := http.Get(url)
	if err != nil {
		panic(err)
	}
	defer res.Body.Close()
	body, err := ioutil.ReadAll(res.Body)
	if err != nil {
		panic(err)
	}
	buf := bytes.NewBuffer(body)
	html := buf.String()
	//fmt.Println(html)
	return html
}

// DNA codon table
func makeDNAcodon() map[string]string {
	codon := make(map[string]string) // DNA codon table
	lines := freadLines("DNAcodon.txt")
	for _, s := range lines {
		line := strings.Split(s, ",")
		codon[line[0]] = line[1]
	}
	return codon
}

func DNAtoProtein(DNA string) string {
	codon := makeDNAcodon()

	prot := ""
	start := strings.Index(DNA, "ATG")
	stop := false
	if start >= 0 {
		for i := start; i < len(DNA)-2; i += 3 {
			v, ok := codon[DNA[i:i+3]]
			if !ok {
				prot = ""
				break
			}
			if v == "Stop" {
				stop = true
				break
			} else {
				prot += v
			}
		}
	}
	if !stop {
		prot = ""
	}
	return prot
}

//stringを貼り合わせてsuper stringを作ることができるという前提条件
//overlapの長さが最も大きいstring同士を貼り合わせて小さいidx要素に移す。用済み要素は空でない末尾要素を移し末尾は空にする。
//有意末尾idxを更新、同じ手順を繰り返す。
func long() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_long.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lt := getIdData(fname)
	lines := []string{}
	//データ行は改行されて続いている
	for i, v := range lt {
		if i%2 == 1 {
			lines = append(lines, v)
		}
	}
	fmt.Println(len(lines))
	// n := len(lines)

	//s1後方とs2前方での一致(前提条件より一方が他方に含まれる場合は除く)
	//ovelapは文字列の長さの半分超
	findCSsub := func(s1, s2 string, lmin int) (int, string) {
		l1 := 0
		ss := ""
		for l := lmin - 1; l > 0; l-- {
			// for l := lmin - 1; l >= lmin/2; l-- {
			flag := 0
			for i := 0; i < l; i++ {
				if s1[len(s1)-l+i] != s2[i] {
					flag = -1
					break
				}
			}
			if flag == 0 {
				ss = s1 + s2[l:]
				l1 = l
				break
			}
		}
		return l1, ss
	}

	findCS := func(s1, s2 string) (int, string) {
		if len(s1) == 0 || len(s2) == 0 {
			return 0, ""
		}

		len1, len2 := len(s1), len(s2)
		lmin := len1
		if len2 < lmin {
			lmin = len2
		}

		//s1後方とs2前方での一致
		l1, cs1 := findCSsub(s1, s2, lmin)
		//s1前方とs2後方での一致
		l2, cs2 := findCSsub(s2, s1, lmin)

		lmax := 0
		cs := ""
		if l1 >= l2 {
			lmax = l1
			cs = cs1
		} else {
			lmax = l2
			cs = cs2
		}
		return lmax, cs
	}

	LCS := func(lines []string, idMax int) int {
		id1, id2 := -1, -1
		lmax := 0
		lcs := ""
		for i := 0; i < idMax; i++ {
			for j := i + 1; j <= idMax; j++ {
				clen, cs := findCS(lines[i], lines[j])
				if clen > lmax {
					lmax = clen
					lcs = cs
					id1, id2 = i, j
				}
			}
		}
		idMax2 := -1
		if lmax > 0 {
			idMax2 = idMax - 1
			lines[id1] = lcs
			for i := id2; i <= idMax2; i++ {
				lines[i] = lines[i+1]
			}
		}
		return idMax2
	}

	idMax := len(lines) - 1
	for idMax > 0 {
		// fmt.Printf("idMax= %d\n", idMax)
		idMax2 := LCS(lines, idMax)
		if idMax2 == -1 {
			fmt.Println("no lcs, quit")
			break
		}
		fmt.Printf("idMax2= %d\n", idMax2)
		idMax = idMax2
	}
	fmt.Println(lines[0])
	ans := []string{lines[0]}
	fwriteLines(output, ans)
}

//bigIntへの変換でint64を求められる
func pmch() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_pmch.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	data := lines[1]
	fmt.Println(data)

	//var calc func(int64) int64   calcで既にint64をオーバーする場合あり
	var calc func(int) *big.Int
	calc = func(n int) *big.Int {
		nB := big.NewInt(int64(n))
		if n == 1 {
			return nB
		}
		return new(big.Int).Mul(nB, calc(n-1)) //n * calc(n-1)
	}

	lt := []string{"A", "U", "C", "G"}
	ct := make([]int, 4)
	for i, v := range lt {
		ct[i] = strings.Count(data, v) //strings.Countの戻り値はint
	}
	fmt.Println(ct)
	if ct[0] != ct[1] || ct[2] != ct[3] {
		panic("dataが不適切です")
	}
	c1 := calc(ct[0])
	c2 := calc(ct[2])
	fmt.Println(c1)
	fmt.Println(c2)
	ans := new(big.Int).Mul(c1, c2)
	fmt.Println(ans)
	ans2 := []string{fmt.Sprint(ans)}
	fwriteLines(output, ans2)
}

func pper() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_pper.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := freadLines(fname)
	lt := strings.Split(lines[0], " ")
	n, _ := strconv.Atoi(lt[0])
	r, _ := strconv.Atoi(lt[1])
	fmt.Println(n, r)
	s := 1
	for i := 0; i < r; i++ {
		s *= n - i
		s = s % 1000000
	}
	fmt.Println(s)
	ans := []string{fmt.Sprint(s)}
	fwriteLines(output, ans)
}

func prob() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_prob.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := freadLines(fname)
	ss := lines[0]
	lt := strings.Split(lines[1], " ")
	pGC := []float64{}
	for i := 0; i < len(lt); i++ {
		x, _ := strconv.ParseFloat(lt[i], 64)
		pGC = append(pGC, x)
	}
	fmt.Println(ss)
	fmt.Println(pGC)
	ans := []string{}
	for j := 0; j < len(pGC); j++ {
		p1 := math.Log10(pGC[j] / 2.0)
		p2 := math.Log10((1.0 - pGC[j]) / 2.0)
		p := 0.0
		for i := 0; i < len(ss); i++ {
			if ss[i] == 'G' || ss[i] == 'C' {
				p += p1
			} else {
				p += p2
			}
		}
		ans = append(ans, fmt.Sprint(p))
	}
	fmt.Println(ans)
	ans2 := strings.Join(ans, " ")
	ans3 := []string{ans2}
	fwriteLines(output, ans3)

}

func sign() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_sign.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := freadLines(fname)
	fmt.Println(lines[0])
	// n, _ := strconv.Atoi(lines[0])
	n := convStoI(lines[0])
	lt := []string{}
	for i := 1; i <= n; i++ {
		lt = append(lt, fmt.Sprint(i))
	}
	fmt.Println(lt)
	ans := permutateStr(lt, len(lt))
	for i := 0; i < n; i++ {
		m := len(ans)
		for j := 0; j < m; j++ {
			ar := make([]string, n)
			copy(ar, ans[j]) //ar := ans[j]では値コピーになっていない
			ar[i] = "-" + ar[i]
			ans = append(ans, ar)
		}
	}
	fmt.Println(ans)
	ans2 := []string{fmt.Sprint(len(ans))}
	for _, v := range ans {
		s := strings.Join(v, " ")
		ans2 = append(ans2, s)
	}
	fwriteLines(output, ans2)
}

func sseq() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_sseq.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)
	ss := lines[1]
	word := lines[3]
	lt := []int{}
	j := 0
	for i := 0; i < len(ss); i++ {
		if j == len(word) {
			break
		}
		if ss[i] == word[j] {
			lt = append(lt, i+1)
			j++
		}
	}
	fmt.Println(lt)
	ans := fmt.Sprint(lt[0])
	for i := 1; i < len(lt); i++ {
		ans += " " + fmt.Sprint(lt[i])
	}
	ans2 := []string{ans}
	fwriteLines(output, ans2)
}

func tran() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_tran.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)
	s1 := lines[1]
	s2 := lines[3]
	dict := map[byte]string{
		'A': "purine",
		'G': "purine",
		'C': "pyrimidine",
		'T': "pyrimidine",
	}
	// for key, value := range dict {
	// 	fmt.Println(key, value)
	// }
	n1 := 0 //transion
	n2 := 0 //transversion
	for i := 0; i < len(s1); i++ {
		if s1[i] == s2[i] {
			continue
		}
		if dict[s1[i]] == dict[s2[i]] {
			n1++
		} else {
			n2++
		}
	}
	ans := float64(n1) / float64(n2)
	fmt.Println(ans)
	ans2 := []string{fmt.Sprint(ans)}
	fwriteLines(output, ans2)

}

func tree() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_tree.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := freadLines(fname)
	fmt.Println(lines)
	n := convStoI(lines[0])
	adj := make([][]int, n)
	for i := 1; i < len(lines); i++ {
		ss := strings.Split(lines[i], " ")
		n1 := convStoI(ss[0]) - 1
		n2 := convStoI(ss[1]) - 1
		adj[n1] = append(adj[n1], n2)
		adj[n2] = append(adj[n2], n1)
	}
	// fmt.Println(adj)

	connected := make([]bool, n)
	p := 0
	connected[0] = true
	for j := 1; j < n; j++ {
		if !isConnected(adj, 0, j) {
			adj[0] = append(adj[0], j)
			adj[j] = append(adj[j], 0)
			p++
		}
		connected[j] = true
	}
	// fmt.Println(connected)
	fmt.Println(p)
	ans := fmt.Sprint(p)
	fwriteLine(output, ans)

}

type LH struct {
	lo int
	hi int
}

func calcMatching(rna string, lo, hi int, dp map[LH]int) int {
	mapping := map[byte]byte{
		'A': 'U',
		'U': 'A',
		'C': 'G',
		'G': 'C',
	}

	characters := hi - lo + 1

	// if there are an odd number of nucleotides,
	// this is an invalid matching.
	if characters%2 == 1 {
		return 0
	}
	// handles tricky edge cases.
	if lo >= hi || lo >= len(rna) || hi < 0 {
		return 1
	}
	// return answer if it is memoized.
	lh := LH{lo, hi}
	_, ok := dp[lh]
	//計算済みの場合は辞書化した値を返す
	if ok {
		return dp[lh]
	} else {
		curr := rna[lo]
		target := mapping[curr]
		acc := 0
		for i := lo + 1; i <= hi; i += 2 {
			// print("lo= %d, i= %d" %(lo, i))
			if rna[i] == target {
				left := calcMatching(rna, lo+1, i-1, dp)
				right := calcMatching(rna, i+1, hi, dp)
				acc += (left * right) % 1000000
			}
		}
		dp[lh] = acc
		return acc
	}
}

func cat() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_cat.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)
	ss := lines[1]
	dp := make(map[LH]int)
	ans := calcMatching(ss, 0, len(ss)-1, dp)
	ans = ans % 1000000
	fmt.Println(ans)
	fwriteLine(output, fmt.Sprint(ans))

}

func cat2() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_cat.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	ss := lines[1]
	fmt.Println(ss)

	count := func(S string) (int, int, int, int) {
		ctA, ctU, ctC, ctG := 0, 0, 0, 0
		for i := 0; i < len(S); i++ {
			switch S[i] {
			case 'A':
				ctA++
			case 'U':
				ctU++
			case 'C':
				ctC++
			case 'G':
				ctG++
			}
		}
		return ctA, ctU, ctC, ctG
	}

	is_complement := func(a, b byte) bool {
		bl := false
		switch a {
		case 'A':
			if b == 'U' {
				bl = true
			}
		case 'U':
			if b == 'A' {
				bl = true
			}
		case 'C':
			if b == 'G' {
				bl = true
			}
		case 'G':
			if b == 'C' {
				bl = true
			}
		}
		return bl
	}

	memo := make(map[string]int) //{} dictionary to memoize solutions in
	memo[""] = 1                 //empty string is already a perfect matching
	var find_noncrossmatch func(string) int
	find_noncrossmatch = func(S string) int {
		if v, ok := memo[S]; ok {
			return v // returns memoized solution if one exists
		}
		ctA, ctU, ctC, ctG := count(S)
		if ctA != ctU || ctC != ctG {
			memo[S] = 0
			return 0 // no perfect matching exists if the number of complementary bases is not equal
		}
		result := 0
		for i := 1; i < len(S); i++ {
			if is_complement(S[0], S[i]) { // look for all posible edges from the first base
				// iteratively compute how many ways can we make a matching with one edge fixed
				result += find_noncrossmatch(S[1:i]) * find_noncrossmatch(S[i+1:])
			}
		}
		result = result % 1000000
		memo[S] = result // memoize solution
		return result
	}

	ans := find_noncrossmatch(ss)
	fmt.Println(ans)
	fwriteLine(output, fmt.Sprint(ans))
}

func main() {
	time1 := time.Now()
	//long()
	//pmch()
	//pper()
	//prob()
	//sign()
	//sseq()
	//tran()
	//tree()
	cat()
	//cat2()

	fmt.Println("elapsed ", time.Now().Sub(time1))
	fmt.Println("End")
}
