import ABC3.Meta.Claim
import ABC3.Interface.GenEll.GaloisRep
import ABC3.Found.GenEll.Lemma31

/-!
# [GenEll] Theorem 3.8 —— Full Special Linear Galois Actions(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.19。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★このスケルトンの位置づけ

`[IUTchIV]` が直接引く [GenEll] の項目は **`Theorem 2.1`(10 回)・`Theorem 3.8`(2 回)・
`Corollary 4.4`(2 回)** など。うち `Corollary 4.3` / `Corollary 4.4`(§4)は
**本定理の上にしか立たない**——すなわち **§4 の 5 件は本定理を経由して
l 捩れ Galois 表現に律速される**。

★**statement を書くだけで語彙が足りない**ので、
`Interface/GenEll/GaloisRep.lean` の `TorsionGaloisRepData` が下敷きになる。
各フィールドは原文 p.19 に**実際に現れる語**と 1 対 1 に対応させた。

## ★★`sorry` の意味(誤読を防ぐ)

本ファイルの `theorem_3_8` は `sorry` である。**これは正しい状態**である——
`Skeleton/` は statement 専用トラックで、証明を付けない場所だからである。

★ただし **`sorry` を消すことを目的にしてはならない**。
`Heights.lean` で実演したとおり、内容を `Interface` の仮説へ移せば `sorry` は消えるが、
それは `tools/check.mjs` 冒頭 B5 が名指しする穴である。
**本当の進捗は `Interface` が `waiting` でなくなること**であり、
`node tools/check.mjs` の「Interface 実装待ち」がそれを見せる。

## ★原文の証明のうち、群論は 1 点だけ(実測)

原文 p.20 の証明を通読した結果、**群論に頼るのは `Lemma 3.1, (iv)` ただ 1 つ**である。
残りはすべて算術側(Faltings 高さ・Tate 曲線・半安定還元・`Lemma 3.7`)。
★`Lemma 3.1` の (i)(ii)(iii) は `Found/GenEll/Lemma31.lean` に**実装済み**(`sorry` 無し)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll

/-- **[GenEll] Theorem 3.8**(Full Special Linear Galois Actions)。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

★**条件 (a) の指数 `ϵ` は `.txt` から読んではならない**——
`pdftotext -layout` は `C · d^ϵ` を `C · d ` と出し、**条件が別物になる**
(2026-08-16 実測。`1_Structured/…/theorem-3-8.html` に記録)。
本 statement は **260 dpi 目視**から写した。

★`23040 = |GL₂(𝔽₃) × GL₂(𝔽₅)|` は原文 p.20 の証明に現れる定数
`(3²−1)(3²−3)(5²−1)(5²−5)` である(目視確認)。 -/
theorem theorem_3_8 (D : TorsionGaloisRepData)
    (KV : Set D.EllClass) (_hKV : D.CompactlyBounded KV)
    (ε : ℝ) (_hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.cls E ∉ Exc →
        ((23040 * 100 * (D.degOfDefinition E : ℝ)
              * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ D.HasPotMultRed E)
          ∨ (D.cls E ∈ KV ∧ D.PrimeToLocalHeights E l ∧ Nat.Coprime l 30)) →
        D.ImageContainsSL2 E l := sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_3_8.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8",
    sectionId := "genell-thm-3-8" }

/-- ★この定理の証明が要求するもの(原文 p.20 を通読して数えた)。

★**群論は `Lemma 3.1, (iv)` ただ 1 つ**で、そのうち (i)(ii)(iii) は実装済みである。
残りはすべて算術側で、そこが `Interface` に載っている。 -/
def theorem_3_8.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Lemma 3.1, (iv)(SL₂(ℤ_l) の持ち上げ)——★原文の証明が使う群論はこれだけ。(i)(ii)(iii) は Found/GenEll/Lemma31.lean に実装済み" 14,
    .otherPaper "[GenEll]" "Lemma 3.7(局所高さと l-巡回部分群スキーム)" 18,
    .otherPaper "[GenEll]" "Proposition 3.4(Faltings 高さによる例外集合の有限性)" 17,
    .otherPaper "[GenEll]" "Lemma 3.2 の直前の局所理論(Tate 曲線)" 15,
    .citation "[FC]" "Degenerations of Abelian Varieties(半安定還元)"
      (.absent "mathlib は `EllipticCurve/Reduction.lean` を持つが Tate 曲線・半安定還元の理論は無い(2026-08-16、ディレクトリ全宣言を確認)") 19,
    .implicitStep
      "原文は 3・5 捩れ点を有理化する次数 23040 の Galois 拡大へ移る段を『there exists a Galois extension L′ of L of degree that divides d₀』と 1 文で済ませている" 20,
    .implicitStep
      "★statement の語彙(Faltings 高さ・Galois 表現・compactly bounded・Galois-finite)を Interface/GenEll/GaloisRep.lean に posit した。**我々は作っていない**" 19 ]

end ABC3.Skeleton.GenEll
