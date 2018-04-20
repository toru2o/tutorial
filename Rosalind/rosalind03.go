/*

func isConnected(adjacent [][]int, start, goal int) bool
func convStoI(s string) int
func freadLines(fname string) []string
func getIdData(fname string) []string
func fwriteLine(fname string, line string)
func fwriteLines(fname string, lines []string)
func MapString(f func(string) string, vs []string) []string
func combi(ss []int, m int) [][]int
func combinate(ss []string, m int) [][]string
// 文字列（スライスに格納）から重複を許してm個を取り出して並べる順列を出力
func permutateDup(ss []string, m int) [][]string
// 文字列（スライスに格納）からm個を取り出して並べる順列を出力
func permutateStr(ss []string, m int) [][]string
// 整数列の順列をスライスで出力
func permutation(n, m int) [][]int
func getBody(url string) string
// DNA codon table
func makeDNAcodon() map[string]string
func DNAtoProtein(DNA string) string

*/
package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
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
func permutateDup(ss []string, m int) [][]string {
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

//stringがreverse complementも含め２回以上あるかどうかチェック
func corr() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_corr.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)
	n := len(lines) / 2

	compl := map[byte]string{
		'A': "T",
		'T': "A",
		'C': "G",
		'G': "C",
	}

	read1 := make([]string, n)
	read2 := make([]string, n)
	m := len(lines[1])
	for i := 0; i < n; i++ {
		read1[i] = lines[2*i+1]
		for j := 0; j < m; j++ {
			read2[i] += compl[read1[i][m-1-j]]
		}
		// fmt.Println(read2[i])
	}
	correct := make([]int, n)
	for i := 0; i < n; i++ {
		if correct[i] >= 1 {
			continue
		}
		for j := i + 1; j < n; j++ {
			if correct[j] >= 1 {
				continue
			}
			if read1[i] == read1[j] || read1[i] == read2[j] {
				correct[i] += 1
				correct[j] += 1
			}
		}
	}
	fmt.Println(correct)

	ans := []string{}
	for i := 0; i < n; i++ {
		if correct[i] > 0 {
			continue
		}
		for j := 0; j < n; j++ {
			if j == i || correct[j] == 0 {
				continue
			}
			ct := 0
			for k := 0; k < m; k++ {
				if ct > 1 {
					break
				}
				if read1[i][k] != read1[j][k] {
					ct++
				}
			}
			if ct == 1 {
				ans = append(ans, read1[i]+"->"+read1[j])
				break
			}

			ct = 0
			for k := 0; k < m; k++ {
				if ct > 1 {
					break
				}
				if read1[i][k] != read2[j][k] {
					ct++
				}
			}
			if ct == 1 {
				ans = append(ans, read1[i]+"->"+read2[j])
				break
			}
		}
	}
	fmt.Println(ans)
	fwriteLines(output, ans)
}

func inod() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_inod.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := freadLines(fname)
	fmt.Println(lines)
	N := convStoI(lines[0])
	//2枚のleafを1個のinternal nodeで繋ぐ
	// n := N / 2 //n個のinternal node(枝1本)を追加。追加したnode同士を繋いでよいが最終的には連結しなければならない
	// r := N % 2 //残ったleaf数(0 or 1)
	//結局leaf相当枚数 n+r を繋いでゆく作業(最終的にはこれらleaf相当は連結されなければならない
	// n+r=2 --> 2つの枝を繋いで終わり(追加nodeなし)
	// n+r=3 --> 2つの枝を繋ぐと(nodeを1個追加) --> n+r=2
	// n+r=4 --> 枝を2つずつ繋いで(nodeを2個追加) --> n+r=2
	node := N / 2
	n := N/2 + (N % 2)
	for n > 2 {
		node += n / 2
		n = n/2 + (n % 2)
		fmt.Println(n)
	}
	if n == 2 {
		fmt.Println(node)
		fwriteLine(output, fmt.Sprint(node))
	} else {
		fmt.Println("no answer")
		fmt.Println(n, node)
	}

}

func kmer() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_kmer.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)

	ss := []string{"A", "C", "G", "T"}
	lt := permutateDup(ss, 4)
	// fmt.Println(lt)

	kmer := make(map[string]int)
	for i := 0; i < len(lt); i++ {
		s := strings.Join(lt[i], "")
		kmer[s] = 0
	}
	// fmt.Println(kmer)

	for i := 0; i < len(lines[1])-3; i++ {
		// fmt.Println(lines[1][i : i+4])
		s := lines[1][i : i+4]
		kmer[s] += 1
	}
	ans := ""
	for _, v := range lt {
		s := strings.Join(v, "")
		ans += fmt.Sprint(kmer[s]) + " "
	}
	fmt.Println(ans)
	fwriteLine(output, ans)

}

func kmp0() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_kmp.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	fmt.Println(lines)
	ss := lines[1]
	P := make([]int, len(ss))

	for i := 1; i < len(ss); i++ {
		for j := i; j > 0; j-- {
			// fmt.Printf("%s, %s\n", ss[i-j+1:i+1], ss[:j])
			if ss[i-j+1:i+1] == ss[:j] {
				P[i] = j
				break
			}
		}
	}
	// fmt.Println(P)
	ans := fmt.Sprint(P[0])
	for i := 1; i < len(P); i++ {
		ans += " " + fmt.Sprint(P[i])
	}
	fmt.Println(ans)
	fwriteLine(output, ans)
}

func kmp() {
	fname := "C:/Users/OGURA/go/practice/data/rosalind_kmp.txt"
	output := "C:/Users/OGURA/go/practice/output.txt"
	lines := getIdData(fname)
	// fmt.Println(lines)
	ss := lines[1]
	T := make([]int, len(ss))

	for i := 1; i < len(ss); i++ {
		for j := 0; i+j < len(ss); j++ {
			if ss[i+j] != ss[j] {
				T[i] = j
				break
			}
		}

	}
	// fmt.Println(T)

	P := make([]int, len(ss))
	for i := 1; i < len(T); i++ {
		if T[i] == 0 {
			continue
		}

		for j := T[i]; j > 0; j-- {
			if T[i] > P[i+j-1] {
				P[i+j-1] = T[i]
			}
			T[i]--
		}
	}
	// fmt.Print(P)

	ans := fmt.Sprint(P[0])
	for i := 1; i < len(P); i++ {
		ans += " " + fmt.Sprint(P[i])
	}
	fmt.Println(ans)
	fwriteLine(output, ans)
}

func main() {
	time1 := time.Now()
	//corr()
	//inod()
	//kmer()
	kmp()

	fmt.Println("elapsed ", time.Now().Sub(time1))
	fmt.Println("End")
}
