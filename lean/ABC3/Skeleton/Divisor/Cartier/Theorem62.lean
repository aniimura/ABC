/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Found.Divisor.CartierMonoid
import ABC3.Skeleton.Divisor.Cartier.Example61

/-!
# Cartier —— `[FrdI] Theorem 6.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta
universe u
variable (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]

/-! ## ★4. Cartier 因子の引き戻し(鎖 `cartier` の `cartier-pullback`)

★★**`Example 6.1` の関手性の本体**である。

## ★★在庫の所在と、いま埋まらない理由(2026-09-06)

★実物は `Found/Divisor/SchemeCartierPull.lean` に sorry ゼロで在る ——
`cartierPullback` / `pullCoeff_add` / `isCartierDiv_cartierPullback` / `pullCoeff_nonneg`。
本節の 4 宣言はこの 4 本にそのまま対応する。

★★**それでも配線できない。** `Found` 側が要求して本節の statement に無いものは 4 つ:

| 足りないもの | どこで要るか |
|---|---|
| `[IsDominant ψ]` | `ffMap ψ : K(X) → K(Y)`(関数体の引き戻し)を作るのに要る |
| `hnormX` / `hnormY` | `ord` が定義できること(正規性が無ければ茎は DVR でない) |
| `hdim : ∀ w, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1` | 局所方程式の取り替えに依らないこと |

★逆に **`[CompactSpace Y]` は足さなくてよい** —— `[AlgebraicGeometry.IsNoetherian Y]` が
局所 Noether ＋ 準コンパクトなので、そこから instance で出る(実測確認済み)。
★`ResearchPaper/decisions-pending.md` の **D8** はここを
「`[IsDominant ψ]`・`hdim`・`[CompactSpace Y]` を足す」と書いていたが、
**`[CompactSpace Y]` は不要**で、代わりに **`hnormX` / `hnormY` が要る**。

★4 つを足したうえでなら、4 宣言は `Found` の 4 本で閉じる(実測済み。
`pullbackCartier` ← `cartierPullback`、`_add` ← `pullCoeff_add`、
`isCartierDiv_` ← `isCartierDiv_cartierPullback`、`_nonneg` ← `pullCoeff_nonneg`)。
★ただし `IsCartierDiv` 自体も `Found` の側(`hnorm` を引数に取る版)へ
移す必要がある —— 本節が使っている `Example61.lean` の `IsCartierDiv` は
`ordAtDiv`(本体が `sorry`)で書かれているからである。 -/

/-- ★★**支配射に沿った Cartier 因子の引き戻し**。

★`Φ(L) → Φ(M)`(`L ⊆ M`)と `Theorem 6.2, (i)` の `Φ₁ → Φ₂|𝒟₁` は、
どちらもこれである。 -/
noncomputable def pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (_ψ : Y ⟶ X)
    (_D : WeilDiv X) (_hD : IsCartierDiv X _D) : WeilDiv Y :=
  sorry

/-- ★引き戻しは加法的。 -/
theorem pullbackCartier_add {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D E : WeilDiv X} (hD : IsCartierDiv X D) (hE : IsCartierDiv X E) :
    pullbackCartier X ψ (D + E) (isCartierDiv_add X hD hE)
      = pullbackCartier X ψ D hD + pullbackCartier X ψ E hE := by
  sorry

/-- ★引き戻しは Cartier 性を保つ。 -/
theorem isCartierDiv_pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D : WeilDiv X} (hD : IsCartierDiv X D) :
    IsCartierDiv Y (pullbackCartier X ψ D hD) := by
  sorry

/-- ★引き戻しは有効性を保つ(`Φ` を `Φ` へ移すために要る)。 -/
theorem pullbackCartier_nonneg {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D : WeilDiv X} (hD : IsCartierDiv X D) (_hpos : 0 ≤ D) :
    0 ≤ pullbackCartier X ψ D hD := by
  sorry

def pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — pulling back divisors",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_add.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しの加法性",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_add.needs : List ProofObligation :=
  [ .derivation "局所的に `f ↦ f ∘ ψ` を取るだけなので、積が和に移る" 110 ]

def isCartierDiv_pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは Cartier 性を保つ",
    sectionId := "frdi-thm-6-2" }

def isCartierDiv_pullbackCartier.needs : List ProofObligation :=
  [ .derivation "局所主因子の引き戻しは局所主因子" 110,
    .implicitStep
      "★原文は「by pulling back [Cartier] divisors and rational functions via ψ, we obtain compatible natural transformations」で畳む" 110 ]

def pullbackCartier_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは有効性を保つ",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_nonneg.needs : List ProofObligation :=
  [ .derivation "支配射なら局所方程式が 0 にならないので、位数は非負のまま" 110 ]

end ABC3.Skeleton.Divisor
