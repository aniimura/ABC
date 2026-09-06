/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Found.Divisor.CartierMonoid
import ABC3.Found.Divisor.SchemeCartierPull
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

## ★★★2026-09-06(判断 D18 採用): `Found` の実物へ配線した

★実物は `Found/Divisor/SchemeCartierPull.lean` に sorry ゼロで在る ——
`cartierPullback` / `pullCoeff_add` / `isCartierDiv_cartierPullback` / `pullCoeff_nonneg`。
本節の 4 宣言はこの 4 本にそのまま対応する。

★★**なぜ待てなかったか(2026-09-06)** —— 2026-09-06 の D8 で `IsCartierDiv` が
`∃ hnorm : IsNormalScheme X, …` の形になったため、旧 statement の

```
IsCartierDiv Y (pullbackCartier X ψ D hD)
```

は「`Y` が正規である」ことまで主張することになった。しかし旧 statement の `Y` は
`[IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y]` しか持たないので、
**非正規な整 Noether スキーム(結節 3 次曲線など)を取れば偽**である。
★**要求を結論に入れることは、要求を消すことではなく義務に変えることである。**
保留を続けると偽の statement が木に残るので、判断 D18 は「(A) で直す」を採った。

## ★★足した仮定は 3 つ。うち逸脱は 1 つだけ

| 足した仮定 | 位置づけ |
|---|---|
| `[IsDominant ψ]` | ★**原典が明示している**(原文 p.110「dominant morphism of schemes」)。**逸脱ではない** |
| `hnormY : IsNormalScheme Y`(源側) | ★原典は `V_i` の**両方**に正規性を常置する("proper normal variety")。**逸脱ではない** |
| `hdim : ∀ w, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1` | ★★**原典に無い。逸脱**(下の逸脱記録を見よ) |

★**的側の正規性 `hnormX` は足していない** —— 仮定 `hD : IsCartierDiv X D` から
`isNormalScheme_of_isCartierDiv` で取り出せる(D8 が正規性を `∃` で述語の内側に畳んだ帰結)。
★**`[CompactSpace Y]` も足していない** —— `[AlgebraicGeometry.IsNoetherian Y]` は
局所 Noether ＋ 準コンパクトなので instance で出る(実測確認済み)。

## ★★★逸脱の記録(判断 D18、2026-09-06)

**`hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1`
は原典に無い仮定である。**

原文は支配射 `ψ` に沿った引き戻しを「by pulling back [Cartier] divisors and rational
functions via ψ」の一言で置き、余次元の保存には触れない。`Found` 側
(`Found/Divisor/SchemeCartierPull.lean`)は**局所方程式の取り替えに依らないこと**を
示す段でこれを要求する —— `f₁/f₂` が `𝒪_{X,ψ(w)}` の単元であることを言うのに、
`ψ(w)` が余次元 1 なら `ord = 0` から、そうでなければ茎が体であることから、
という 2 分岐を使うためである。
★`Found` 側には以前から入っていたが、**逸脱として記録された形跡が無かった**ので
ここで記録する(D18 の測り直しで見つかった記録の穴の 1 つ)。
★原典の 2 つの用例(正規化射 `V[M] → V[L]`・Frobenius)では自動的に成立するので
臨界路には乗らない。局所 Hartogs で `hdim` を落とす道は**独立ノードとして後送り**する。

## ★橋は 2 本だけ

`Skeleton` の `IsCartierDiv X D`(`∃ hnorm, …`、局所方程式は `Kˣ`)と
`Found` の `IsCartierDiv hnorm D`(局所方程式は `(f, f ≠ 0)`)の間の詰め替えである。
`ordAtDiv = ordPt` は定義なので、実体は `f : Kˣ ↔ (f, f ≠ 0)` の往復だけであり、
`hnorm` の取り替えは `IsNormalScheme X : Prop` の**証明無関係性**で defeq になる。 -/

/-- ★`Skeleton` の `IsCartierDiv` は正規性を `∃` で抱えているので、そこから取り出せる。 -/
theorem isNormalScheme_of_isCartierDiv {D : WeilDiv X} (hD : IsCartierDiv X D) :
    IsNormalScheme X :=
  hD.elim fun h _ => h

/-- ★★**橋(→)** —— `Skeleton` の Cartier 性から `Found` の Cartier 性へ。

★`hnorm` は任意に取ってよい(`IsNormalScheme X` は `Prop` なので証明無関係)。 -/
theorem found_isCartierDiv_of_isCartierDiv {D : WeilDiv X} (hD : IsCartierDiv X D)
    (hnorm : IsNormalScheme X) : ABC3.Found.Divisor.IsCartierDiv hnorm D := by
  obtain ⟨_, h⟩ := hD
  intro x
  obtain ⟨U, hxU, f, hf⟩ := h x
  exact ⟨U, hxU, (f : X.functionField), Units.ne_zero f, fun v hv => hf v hv⟩

/-- ★★**橋(←)** —— `Found` の Cartier 性から `Skeleton` の Cartier 性へ。 -/
theorem isCartierDiv_of_found {D : WeilDiv X} (hnorm : IsNormalScheme X)
    (hD : ABC3.Found.Divisor.IsCartierDiv hnorm D) : IsCartierDiv X D := by
  refine ⟨hnorm, fun x => ?_⟩
  obtain ⟨U, hxU, f, hfne, hf⟩ := hD x
  exact ⟨U, hxU, Units.mk0 f hfne, fun v hv => hf v hv⟩

/-- ★★**支配射に沿った Cartier 因子の引き戻し**。

★`Φ(L) → Φ(M)`(`L ⊆ M`)と `Theorem 6.2, (i)` の `Φ₁ → Φ₂|𝒟₁` は、
どちらもこれである。

★★★2026-09-06(D18)に `Found` へ配線した(`ABC3.Found.Divisor.cartierPullback`)。
足した仮定は `[IsDominant ψ]`・`hnormY`・`hdim` の 3 つ
(前 2 つは原典の仮定の復元、`hdim` は逸脱。冒頭の逸脱記録を見よ)。 -/
noncomputable def pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) [IsDominant ψ]
    (hnormY : IsNormalScheme Y)
    (hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1)
    (D : WeilDiv X) (hD : IsCartierDiv X D) : WeilDiv Y :=
  ABC3.Found.Divisor.cartierPullback ψ hnormY (isNormalScheme_of_isCartierDiv X hD)
    (found_isCartierDiv_of_isCartierDiv X hD (isNormalScheme_of_isCartierDiv X hD)) hdim

/-- ★成分は `Found` の `pullCoeff`(局所方程式の引き戻しの位数)である。 -/
theorem pullbackCartier_apply {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) [IsDominant ψ]
    (hnormY : IsNormalScheme Y)
    (hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1)
    (D : WeilDiv X) (hD : IsCartierDiv X D) (w : PrimeDivisorPt Y) :
    pullbackCartier X ψ hnormY hdim D hD w
      = ABC3.Found.Divisor.pullCoeff ψ hnormY (isNormalScheme_of_isCartierDiv X hD)
          (found_isCartierDiv_of_isCartierDiv X hD (isNormalScheme_of_isCartierDiv X hD)) w :=
  rfl

/-- ★引き戻しは加法的。

★★★2026-09-06(D18)に `Found` へ配線した(`ABC3.Found.Divisor.pullCoeff_add`)。 -/
theorem pullbackCartier_add {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) [IsDominant ψ]
    (hnormY : IsNormalScheme Y)
    (hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1)
    {D E : WeilDiv X} (hD : IsCartierDiv X D) (hE : IsCartierDiv X E) :
    pullbackCartier X ψ hnormY hdim (D + E) (isCartierDiv_add X hD hE)
      = pullbackCartier X ψ hnormY hdim D hD + pullbackCartier X ψ hnormY hdim E hE := by
  refine Finsupp.ext (fun w => ?_)
  rw [Finsupp.add_apply, pullbackCartier_apply, pullbackCartier_apply, pullbackCartier_apply]
  exact ABC3.Found.Divisor.pullCoeff_add ψ hnormY (isNormalScheme_of_isCartierDiv X hD)
    (found_isCartierDiv_of_isCartierDiv X hD _) (found_isCartierDiv_of_isCartierDiv X hE _)
    (found_isCartierDiv_of_isCartierDiv X (isCartierDiv_add X hD hE) _) w (hdim w)

/-- ★引き戻しは Cartier 性を保つ。

★★★2026-09-06(D18)に `Found` へ配線した
(`ABC3.Found.Divisor.isCartierDiv_cartierPullback`)。
★結論の `IsCartierDiv Y (…)` は D8 以後「`Y` の正規性」も主張するので、
`hnormY` が**必須**になっている(これが D18 を待てなくした理由である)。 -/
theorem isCartierDiv_pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) [IsDominant ψ]
    (hnormY : IsNormalScheme Y)
    (hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1)
    {D : WeilDiv X} (hD : IsCartierDiv X D) :
    IsCartierDiv Y (pullbackCartier X ψ hnormY hdim D hD) :=
  isCartierDiv_of_found Y hnormY
    (ABC3.Found.Divisor.isCartierDiv_cartierPullback ψ hnormY
      (isNormalScheme_of_isCartierDiv X hD)
      (found_isCartierDiv_of_isCartierDiv X hD (isNormalScheme_of_isCartierDiv X hD)) hdim)

/-- ★引き戻しは有効性を保つ(`Φ` を `Φ` へ移すために要る)。

★★★2026-09-06(D18)に `Found` へ配線した(`ABC3.Found.Divisor.pullCoeff_nonneg`)。
`Finsupp` の順序(`Finsupp.le_def`)を成分ごとの不等式に開くところだけが橋である。 -/
theorem pullbackCartier_nonneg {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) [IsDominant ψ]
    (hnormY : IsNormalScheme Y)
    (hdim : ∀ w : PrimeDivisorPt Y, ringKrullDim (X.presheaf.stalk (ψ.base w.1)) ≤ 1)
    {D : WeilDiv X} (hD : IsCartierDiv X D) (hpos : 0 ≤ D) :
    0 ≤ pullbackCartier X ψ hnormY hdim D hD := by
  have hDnn : ∀ v : PrimeDivisorPt X, 0 ≤ D v := by
    intro v
    simpa using (Finsupp.le_def.mp hpos) v
  refine Finsupp.le_def.mpr (fun w => ?_)
  rw [Finsupp.coe_zero, Pi.zero_apply, pullbackCartier_apply]
  exact ABC3.Found.Divisor.pullCoeff_nonneg ψ hnormY (isNormalScheme_of_isCartierDiv X hD)
    (found_isCartierDiv_of_isCartierDiv X hD _) w (hdim w) hDnn

def isNormalScheme_of_isCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Cartier 性から正規性を取り出す(Skeleton ↔ Found の橋)",
    sectionId := "frdi-thm-6-2" }

def isNormalScheme_of_isCartierDiv.needs : List ProofObligation :=
  [ .derivation "`Skeleton` の `IsCartierDiv` は D8 以後 `∃ hnorm : IsNormalScheme X, …` の形なので、`Exists.elim` で取り出せる" 110 ]

def found_isCartierDiv_of_isCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Skeleton の Cartier 性から Found の Cartier 性へ(橋)",
    sectionId := "frdi-thm-6-2" }

def found_isCartierDiv_of_isCartierDiv.needs : List ProofObligation :=
  [ .citation "[ABC3]" "IsCartierDiv(Found 側の定義。局所方程式は `(f, f ≠ 0)`)"
      (.inProject "ABC3" "ABC3.Found.Divisor.IsCartierDiv") 110,
    .derivation "`ordAtDiv = ordPt` は定義なので、詰め替えるのは局所方程式 `f : Kˣ ↔ (f, f ≠ 0)` だけ。`hnorm` の取り替えは `IsNormalScheme X : Prop` の証明無関係性で defeq" 110 ]

def isCartierDiv_of_found.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — Found の Cartier 性から Skeleton の Cartier 性へ(橋)",
    sectionId := "frdi-thm-6-2" }

def isCartierDiv_of_found.needs : List ProofObligation :=
  [ .derivation "`Units.mk0` で `f ≠ 0` を `Kˣ` に戻し、`∃ hnorm` の証人に `hnorm` を置く" 110 ]

def pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — pulling back divisors",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_apply.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 引き戻しの成分は局所方程式の引き戻しの位数",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_apply.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cartierPullback_apply(Found 側。`Finsupp.onFinset` の成分は `pullCoeff`)"
      (.inProject "ABC3" "ABC3.Found.Divisor.cartierPullback_apply") 110,
    .derivation "`Finsupp.onFinset` の成分は定義そのものなので `rfl`" 110 ]

/-- ★★2026-09-06(D18)に新設。`pullbackCartier` には `.needs` が無かった。

★原文が「by assumption (a)」と名指ししている**台の条件**をここに写す
(D18 の測り直しで見つかった記録の穴)。 -/
def pullbackCartier.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cartierPullback(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.cartierPullback") 110,
    .citation "[ABC3]" "ordPt_ffMap_congr(引き戻しの係数は局所方程式の取り替えに依らない)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ordPt_ffMap_congr") 110,
    .citation "[ABC3]" "pullCoeff_eq_zero_of_notMem(原文の仮定 (a)——D_{K₁} の引き戻しが D_{K₂} に入ること——から引き戻しの台が D_{L₂} に収まる)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_eq_zero_of_notMem") 110,
    .derivation "支配射なので `ffMap ψ : K(X) → K(Y)` が存在し、局所方程式 `f` の引き戻し `ψ^* f` の位数で係数を定める" 110,
    .implicitStep
      "★原文は「by assumption (a)」の一言で台の条件を引く。本節の statement には台の条件を持ち込まず、Found の Thm62Pull.lean が `hpull` として抱えている(底が動く版)" 110,
    .implicitStep
      ("★逸脱(D18、2026-09-06): `hdim`(ψ の像の茎の Krull 次元 ≤ 1)は原典に無い仮定である。"
       ++ "局所方程式の取り替えに依らないことを示す段で Found 側が要求する。"
       ++ "原典の 2 用例(正規化射・Frobenius)では自動成立するので臨界路には乗らない") 110 ]

def pullbackCartier_add.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しの加法性",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_add.needs : List ProofObligation :=
  [ .derivation "局所的に `f ↦ f ∘ ψ` を取るだけなので、積が和に移る" 110,
    .citation "[ABC3]" "pullCoeff_add(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_add") 110 ]

def isCartierDiv_pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは Cartier 性を保つ",
    sectionId := "frdi-thm-6-2" }

def isCartierDiv_pullbackCartier.needs : List ProofObligation :=
  [ .derivation "局所主因子の引き戻しは局所主因子" 110,
    .citation "[ABC3]" "isCartierDiv_cartierPullback(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isCartierDiv_cartierPullback") 110,
    .implicitStep
      "★原文は「by pulling back [Cartier] divisors and rational functions via ψ, we obtain compatible natural transformations」で畳む" 110 ]

def pullbackCartier_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは有効性を保つ",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_nonneg.needs : List ProofObligation :=
  [ .derivation "支配射なら局所方程式が 0 にならないので、位数は非負のまま" 110,
    .citation "[ABC3]" "pullCoeff_nonneg(Found 側の本体、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullCoeff_nonneg") 110 ]

end ABC3.Skeleton.Divisor
