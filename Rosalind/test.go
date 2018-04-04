package main

import (
	"bufio"
	"fmt"
	// "math/big"
	"os"
)

// const BUFSIZE = 1024

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
		fmt.Printf("%4d行目: %s\n", i, sc.Text())
		lines = append(lines, sc.Text())
	}
	return lines
}

func fwriteLines(fname string, lines []string) {
	fout, err := os.Create(fname)
	if err != nil {
		os.Exit(1)
	}

	for i := 0; i < len(lines); i++ {
		fout.Write(([]byte)(lines[i] + "\r\n"))
	}
	fout.Close()
}

type Listi []int

func (xs Listi) Sum() int {
	s := 0
	for _, v := range xs {
		s += v
	}
	return s
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

type LH struct {
	lo int
	hi int
}

func main() {
	memo := make(map[string]int)
	memo["ab"] = 1
	memo["cd"] = 2
	_, ok := memo["ab"]
	if ok {
		fmt.Println(memo["ab"])
	}
	if _, ok := memo["cd"]; ok {
		fmt.Println(memo["cd"])
	}
	if v, ok := memo["cd"]; ok {
		fmt.Println(v)
	}
	if v, ok := memo["ef"]; ok {
		fmt.Println(v)
	} else {
		fmt.Println("not exists")
	}

}
