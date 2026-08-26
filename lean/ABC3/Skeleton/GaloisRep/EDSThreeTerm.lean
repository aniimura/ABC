import ABC3.Found.GaloisRep.EDSThreeTerm

/-!
# スケルトン —— **`normEDS` の 3 項恒等式**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★(G5) の最下流の葉

(G5) の非退化性 → `deg[n] = n²`(第 196・197)→ `x([n]P) = Φ_n/ΨSq_n`(第 198)
→ **EDS 恒等式**、という縮約の**最後の的**である。

★第 322 ブロック(`Found/GaloisRep/EDSThreeTerm.lean`)で、
`IsEllSequence` は **2 変数の 3 項恒等式**から `ring` で出ることを示した。
★★したがって残るのは本ファイルの `normEDS_isThreeTerm` 1 件だけである。

## ★★★★★既に分かっている場合

| 場合 | 出どころ | 状態 |
|---|---|---|
| `n = 0`・`n = 1` | 自明 | ✅ `normEDS_threeTerm_zero`・`_one` |
| 差 `m − n = 1` | `normEDS_odd` | ✅ `normEDS_threeTerm_diff_one` |
| 差 `m − n = 2` | `normEDS_even` | ✅ `normEDS_threeTerm_diff_two` |
| `m ↔ n` の交換 | 奇関数性 | ✅ `threeTerm_swap` |
| **差 `3` 以上** | —— | ★**これだけが残る** |

## ★★★★★★差 `3` 以上が「局所」には出ない理由

mathlib の 2 本の漸化式は指数 `n` を `n/2` に落とす。★したがって差 `3` の場合

    c·[W(m)W(m−2)³ − W(m−3)W(m−1)³] = W(m+1)W(m−1)W(m−3)² − W(m−2)W(m−4)W(m)²

を隣接指数だけの `ring` で出すことはできない——**倍加の構造を経由する**必要がある。
★★実際、第 202 で `n = 3` の梯子恒等式は `ring` で通ったが、
第 203 相当の `n = 4`(`n → 2n` の段)は 124 秒かけて失敗した(2026-08-20 実測)。

## ★★これは mathlib 自身の TODO である

`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean`:
> TODO: prove that `normEDS` satisfies `IsEllDivSequence`.

★★★★上流に入れるのが本筋である。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep

variable {R : Type} [CommRing R]

/-! ## ★★★★★★`normEDS` は 3 項恒等式を満たす -/

/-- ★★★★★★**`normEDS` は 3 項恒等式を満たす**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 322 で `IsEllSequence` はここから `ring` で出る(`isEllSequence_of_isThreeTerm`)。
★★差 `1`・差 `2` は mathlib の漸化式そのもの、`n = 0, 1` は自明。
★★★残るのは**差 `3` 以上**であり、倍加の構造を経由した強い帰納法が要る。 -/
theorem normEDS_isThreeTerm (b c d : R) : IsThreeTerm (normEDS b c d) := by
  sorry

/-- ★★★★★**`normEDS` は楕円列である**——第 322 の縮約を通す。 -/
theorem normEDS_isEllSequence (b c d : R) : IsEllSequence (normEDS b c d) :=
  isEllSequence_of_isThreeTerm (normEDS_isThreeTerm b c d)

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def normEDS_isEllSequence.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分点多項式列が楕円列であること)",
    sectionId := "genell-thm-3-8" }

def normEDS_isEllSequence.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**第 322 ブロックの縮約を通すだけである**——`isEllSequence_of_isThreeTerm`(`ring` で閉じる)に `normEDS_isThreeTerm` を当てる。★したがって本節点に固有の数学は無く、上の `normEDS_isThreeTerm` に完全に従属する(0 ブロック)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の非退化性——分点多項式列が楕円列であること)" 19 ]

def normEDS_isThreeTerm.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分点多項式列が楕円列であること)",
    sectionId := "genell-thm-3-8" }

def normEDS_isThreeTerm.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**2026-08-26(第 322 ブロック): 的が 3 変数から 2 変数に縮んだ**。`IsEllSequence` の主張 `h(m,n)h(r,0) = h(m,r)h(n,0) − h(n,r)h(m,0)`(`h(x,y) := W(x+y)W(x−y)`)は、`r = 1` の場合が `h(x,y) = det((W(x+1)W(x−1)), (W(x)²))` という **2 次元ベクトルの行列式表示**を与えるので、あとは **2 次元の Plücker 関係**として `ring` で閉じる(`isEllSequence_of_isThreeTerm`)。★Riemann の 4 項恒等式も同じ理由で出る(0 ブロック)" 19,
    .implicitStep
      "★★★★★**2026-08-26: 既知の場合を洗い出した**。`n = 0`・`n = 1` は自明(`normEDS_threeTerm_zero`・`_one`)、差 `m − n = 1` は `normEDS_odd` そのもの(`normEDS_threeTerm_diff_one`)、差 `2` は `normEDS_even` そのもの(`normEDS_threeTerm_diff_two`)、`m ↔ n` の交換は奇関数性(`threeTerm_swap`)。★★**残るのは差 `3` 以上だけ**である(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**残件: 差 `3` 以上**。mathlib の漸化式は指数を半分に落とすので、隣接指数だけの `ring` では出ない——**倍加の構造を経由した強い帰納法**が要る。★偶奇 4 通りの場合分けで `W(2a±2b)` を `W(a±b+k)`(`|k| ≤ 2`)に展開し、帰納法の仮定を半分の指数の対に当てる形になる。★★見積もり 30-80 ブロック" 19,
    .citation "[Ward]" "Memoir on Elliptic Divisibility Sequences(楕円列の 3 項恒等式)"
      (.absent "mathlib は `IsEllSequence`・`normEDS` を持つが、`normEDS` が `IsEllDivSequence` を満たすことは **TODO のまま**(2026-08-26 実測、`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean` の Main statements 節)") 19,
    .implicitStep
      "★★割り切れ性(`IsDivSequence`)の側は mathlib が `complEDS` を用意しており、`W(k)·Wᶜ(k,n) = W(nk)` が示せれば出る。★`n = 2` の場合(`normEDS_mul_complEDS₂`)だけは**証明済み**(2026-08-26 実測)。一般の `n` は未証明で、これも 3 項恒等式に依る見込み" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の非退化性——分点多項式列が楕円列であること)" 19 ]

end ABC3.Skeleton.GaloisRep
