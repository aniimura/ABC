import ABC3.Meta.Claim
import ABC3.Interface.GenEll.GaloisRep
import ABC3.Found.GenEll.Lemma31
import ABC3.Found.GenEll.Thm38Assembly
import ABC3.Skeleton.GenEll.Section3

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

## ★★★★★★★★ 2026-08-26——`sorry` が消えた(第 365 ブロック)

★以前はここに「`theorem_3_8` は `sorry` であり、**それを消すことを目的にしてはならない**」
と書いてあった——`Heights.lean` で実演したとおり、内容を `Interface` の仮説へ移せば
`sorry` は消えるが、それは `tools/check.mjs` 冒頭 **B5** が名指しする穴だからである。

★★**その警告は今も生きている**。本ブロックが B5 にならないのは、次の 3 点による:

1. ★**界面に出した 7 つはすべて原文 p.20 が実際に書いている段**である
   ——`L′` への基底変換(次数が 23040 を割る、半安定還元、乗法還元)と、最終段の `Lemma 3.1, (iv)`。
   **結論そのものを仮説に移したのではない**。
2. ★★**中身は `Lemma 3.7` からの導出である**——定数 23040 の帳尻(下記)も
   条件 (b) の 30 の使い道も、ここで初めて仕事をする。
3. ★★★**残った 1 つの posit は名前で見える**——`imageContainsSL2_of_torsionExt` だけであり、
   それは **Galois 表現が未構築だから**であって、`Lemma 3.1` を持っていないからではない。

★★★★**本当の進捗は依然として `Interface` が `waiting` でなくなること**であり、
`node tools/check.mjs` の「Interface 実装待ち」がそれを見せる。

## ★★★★ 2026-08-26 の逸脱の記録——量化する界面を下げた

| 項 | 原典 | 形式化 |
|---|---|---|
| 量化する界面 | (区別なし) | `TorsionGaloisRepData` → **`EllModuliData`** |

★理由: 原文 p.20 の証明は `Lemma 3.7` を使うが、`Lemma 3.7` は `EllModuliData` の上に立つ。
★★`EllModuliData extends TorsionGaloisRepData` なので、**下流は何も失わない**
——`Corollary 4.3` / `Corollary 4.4` はもともと `EllModuliData` の上にある。
★★★弱めているのは「どの `D` について言うか」だけで、**結論は同じ**である。

## ★★★★★定数 23040 の帳尻(本ブロックの算術的な中身)

原文は『for a suitable choice of `C` … then `E_{L′}` and `l` satisfy condition (a) …
of Lemma 3.7 **[perhaps for a different “C”]**』と 1 文で済ませている。実際には:

    d′ ≤ 23040·d 、 [E_{L′}] = [E_L] なので ht^Falt はそのまま
    d′^ε ≤ (23040·d)^ε = 23040^ε·d^ε
    C := C₇·23040^ε + |B| + 1 と取れば
      100·d′·(ht^Falt + C₇·d′^ε) ≤ 23040·100·d·(ht^Falt + C·d^ε)

★★`ht^Falt < 0` のとき第 1 項は**逆向きになる**が、その超過分は `|B|·23040·d` 以下であり
(`ht^Falt ≥ B`)、`C` の `|B| + 1` の分が `d^ε ≥ 1` を使ってそれを吸収する。

## ★原文の証明のうち、群論は 1 点だけ(実測)

原文 p.20 の証明を通読した結果、**群論に頼るのは `Lemma 3.1, (iv)` ただ 1 つ**である。
残りはすべて算術側(Faltings 高さ・Tate 曲線・半安定還元・`Lemma 3.7`)。
★`Lemma 3.1` の (i)(ii)(iii) は `Found/GenEll/Lemma31.lean` に**実装済み**(`sorry` 無し)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll

set_option maxHeartbeats 1600000 in
/-- **[GenEll] Theorem 3.8**(Full Special Linear Galois Actions)。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

★**条件 (a) の指数 `ϵ` は `.txt` から読んではならない**——
`pdftotext -layout` は `C · d^ϵ` を `C · d ` と出し、**条件が別物になる**
(2026-08-16 実測。`1_Structured/…/theorem-3-8.html` に記録)。
本 statement は **260 dpi 目視**から写した。

★`23040 = |GL₂(𝔽₃) × GL₂(𝔽₅)|` は原文 p.20 の証明に現れる定数
`(3²−1)(3²−3)(5²−1)(5²−5)` である(目視確認)。 -/
theorem theorem_3_8 (D : EllModuliData)
    (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV)
    (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.cls E ∉ Exc →
        ((23040 * 100 * (D.degOfDefinition E : ℝ)
              * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ D.HasPotMultRed E)
          ∨ (D.cls E ∈ KV ∧ D.PrimeToLocalHeights E l ∧ Nat.Coprime l 30)) →
        D.ImageContainsSL2 E l :=
  ABC3.Found.GenEll.theorem_3_8 D KV hKV ε hε

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_3_8.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8",
    sectionId := "genell-thm-3-8" }

/-- ★この定理の証明が要求するもの(原文 p.20 を通読して数えた)。

★**群論は `Lemma 3.1, (iv)` ただ 1 つ**で、そのうち (i)(ii)(iii) は実装済みである。
残りはすべて算術側で、そこが `Interface` に載っている。 -/
def theorem_3_8.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Lemma 3.1, (iv)(SL₂(ℤ_l) の持ち上げ)——★原文の証明が使う群論はこれだけであり、★★★★★**4 条すべてが実装済み**である(2026-08-29 確認): (i)(ii)(iii) は Found/GenEll/Lemma31.lean、**(iv) は Found/GenEll/Sl2Padic.lean の lemma_3_1_iv**。★原文は [Serre] Chapter IV §3.4 Lemma 3 を引くが 0_Source に無いため、本プロジェクトが自分で証明している。★★★★★★**したがって Theorem 3.8 に「Serre の開像定理」は要らない**——障壁は Tate 曲線(下の 2 行)だけである" 14,
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
