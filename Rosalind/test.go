package main

import (
	"bufio"
	"fmt"
	"math/big"
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

func main() {
	// u, err := url.Parse("http://www.uniprot.org/uniprot/B5ZC00.fasta")
	// if err != nil {
	// 	log.Fatal(err)
	// }
	// fmt.Printf("URL: %s\n", u.String())
	// fmt.Printf("Scheme: %s\n", u.Scheme)
	// fmt.Printf("Opaque: %s\n", u.Opaque)
	// fmt.Printf("User: %s\n", u.User)
	// fmt.Printf("Host: %s\n", u.Host)
	// fmt.Printf("Path: %s\n", u.Path)
	// fmt.Printf("RawPath: %s\n", u.RawPath)
	// fmt.Printf("RawQuery: %s\n", u.RawQuery)
	// fmt.Printf("Fragment: %s\n", u.Fragment)

	// for key, values := range u.Query() {
	// 	fmt.Printf("Query Key: %s\n", key)
	// 	for i, v := range values {
	// 		fmt.Printf("Query Value[%d]: %s\n", i, v)
	// 	}
	// }

	// url := "http://www.uniprot.org/uniprot/B5ZC00.fasta"
	// res, err := http.Get(url)
	// if err != nil {
	// 	panic(err)
	// }
	// defer res.Body.Close()
	// body, err := ioutil.ReadAll(res.Body)
	// if err != nil {
	// 	panic(err)
	// }
	// buf := bytes.NewBuffer(body)
	// html := buf.String()
	// fmt.Println(html)

	x := big.NewInt(2)
	fmt.Println(x)
	fmt.Println("Git test02")

	fmt.Println("Git test03")

}
