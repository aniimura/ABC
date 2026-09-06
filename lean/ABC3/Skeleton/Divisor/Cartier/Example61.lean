/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Found.Divisor.CartierMonoid
import ABC3.Found.Divisor.SchemeCartier

/-!
# Cartier —— `[FrdI] Example 6.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta
universe u
variable (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]

/-! ## ★1. Cartier 因子(鎖 `cartier` の `cartier-def`) -/

/-- ★**Cartier 因子** —— 局所的に有理函数の因子であるもの。

★各点のまわりで、ある `f ∈ K(X)^×` があって `D` がその近傍で `div(f)` に一致する。

## ★★逸脱の記録(2026-09-06、D8 採用)

**2026-09-06(D8 採用): 原典の proper normal variety の正規性を復元した。
Skeleton が落としていた仮定を戻すもので、原典からの逸脱ではない。
`[IsLocallyNoetherian X]` は `ordPt` の要求。`[CompactSpace X]` は `IsNoetherian` から出るので
足していない。**

★★**正規性を引数ではなく `∃` で内側に持たせた**。理由は配線上の制約である ——
`ordAtDiv` が `hnorm` を取るようになったので `IsCartierDiv` も正規性を要するが、
`IsCartierDiv` を `(hnorm) (D)` の 2 引数にすると
`Skeleton/Divisor/Cartier/Theorem62.lean`(**D18 で保留中**なので触れない)の
`IsCartierDiv X D` / `isCartierDiv_add X hD hE` が壊れる。
★`∃ hnorm` にすると引数の数が変わらないので `Theorem62.lean` は無傷のまま通る。

★★**向きは「強める」側である** —— 非正規な `X` では
`IsCartierDiv X D` は**どの `D` についても偽**になる(空虚に真にはならない)。
`∀ hnorm, …` と書けば非正規の枝が空虚に真になり `cartierSubgroup = ⊤` という
**新しい退化**を作るので、そちらは採らなかった。
★原典 `Example 6.1` は `V` を proper normal variety と置いた**上で** Cartier 因子を語るので、
正規性を述語に畳み込むのは原典の常置仮定を写したものである。

★★**副作用の記録**: `Theorem62.lean::isCartierDiv_pullbackCartier` の結論
`IsCartierDiv Y (…)` は、この変更により `Y` の正規性も主張することになった。
D18 が「`hnormY` を足す必要がある」と測定済みの箇所と**同じもの**であり、
当該 4 宣言は現在 `sorry` のままなので新たな偽の証明は生じない。 -/
def IsCartierDiv (D : WeilDiv X) : Prop :=
  ∃ hnorm : IsNormalScheme X, ∀ x : X, ∃ (U : X.Opens) (_ : x ∈ U) (f : (X.functionField)ˣ),
    ∀ y : PrimeDivisorPt X, y.1 ∈ U → D y = ordAtDiv X hnorm y (f : X.functionField)

/-- ★**`Q`-Cartier 因子** —— ある正の倍数が Cartier。 -/
def IsQCartierDiv (D : WeilDiv X) : Prop :=
  ∃ n : ℕ, 0 < n ∧ IsCartierDiv X (n • D)

/-! ## ★2. Cartier 因子は部分群をなす(鎖 `cartier` の `cartier-subgroup`)

★★**ここが単系論への入口である。**`Found/Divisor/CartierMonoid.lean` は
`Γ : AddSubgroup (S →₀ ℤ)` しか要求しないので、この 3 条が出れば
`Example 6.1` の単系論(divisorial・`Φ^gp ≃ Γ`・`Prime ≃ D_L`・perf-factorial)は
**すべて自動で付く**。

## ★★★2026-09-06(D8 採用): `Found` の実物へ配線した

★実物は `Found/Divisor/SchemeCartier.lean` の `cartierSubgroup` に sorry ゼロで在り、
`zero_mem'` / `add_mem'` / `neg_mem'` の 3 条が `ordPt_one` / `ordPt_mul` / `ordPt_inv`
だけで証明されている。本節の 3 条はその 3 行をそのまま写したものである。

★★**`hnorm` を要るのは `isCartierDiv_zero` だけ**である ——
`IsCartierDiv` が正規性を `∃` で抱えているので、`add` / `neg` は
**仮定 `hD` から正規性を取り出せる**。したがって 2 本には追加の仮定を置いていない
(置くと `Theorem62.lean` の `isCartierDiv_add X hD hE` が壊れる)。 -/

/-- ★`0` は Cartier。

★★★2026-09-06 に埋まった。`Found` の `cartierSubgroup.zero_mem'` と同じ 1 行
(`U = ⊤`・`f = 1`・`ordPt_one`)。★`hnorm` を要求するのはこの 1 本だけである。 -/
theorem isCartierDiv_zero (hnorm : IsNormalScheme X) : IsCartierDiv X 0 :=
  ⟨hnorm, fun _ => ⟨⊤, trivial, 1, fun y _ => by
    simp [ordAtDiv, ABC3.Found.Divisor.ordPt_one]⟩⟩

/-- ★Cartier 因子の和は Cartier。

★★★2026-09-06 に埋まった。`Found` の `cartierSubgroup.add_mem'` と同じ
(`U₁ ⊓ U₂` の上で `f₁ * f₂`、`ordAtDiv_mul`)。 -/
theorem isCartierDiv_add {D E : WeilDiv X} (hD : IsCartierDiv X D) (hE : IsCartierDiv X E) :
    IsCartierDiv X (D + E) := by
  obtain ⟨hnorm, hD⟩ := hD
  obtain ⟨_, hE⟩ := hE
  refine ⟨hnorm, fun x => ?_⟩
  obtain ⟨U₁, hx₁, f₁, h₁⟩ := hD x
  obtain ⟨U₂, hx₂, f₂, h₂⟩ := hE x
  refine ⟨U₁ ⊓ U₂, ⟨hx₁, hx₂⟩, f₁ * f₂, fun y hy => ?_⟩
  rw [Finsupp.add_apply, h₁ y hy.1, h₂ y hy.2]
  exact (ordAtDiv_mul X hnorm y _ _ (Units.ne_zero f₁) (Units.ne_zero f₂)).symm

/-- ★Cartier 因子の逆元は Cartier。

★★★2026-09-06 に埋まった。`Found` の `cartierSubgroup.neg_mem'` と同じ
(同じ `U` の上で `f⁻¹`、`ordPt_inv`)。 -/
theorem isCartierDiv_neg {D : WeilDiv X} (hD : IsCartierDiv X D) :
    IsCartierDiv X (-D) := by
  obtain ⟨hnorm, hD⟩ := hD
  refine ⟨hnorm, fun x => ?_⟩
  obtain ⟨U, hx, f, h⟩ := hD x
  refine ⟨U, hx, f⁻¹, fun y hy => ?_⟩
  rw [Finsupp.neg_apply, h y hy]
  simp only [ordAtDiv, Units.val_inv_eq_inv_val]
  rw [ABC3.Found.Divisor.ordPt_inv hnorm y (Units.ne_zero f)]

/-- ★★**Cartier 因子の群** —— `Example 6.1` の `Φ(L)^gp ⊆ ℤ[D_L]` の中身。

★`hnorm` は `zero_mem'` のためだけに要る(`add_mem'` / `neg_mem'` は仮定から取り出せる)。 -/
def cartierSubgroup (hnorm : IsNormalScheme X) : AddSubgroup (WeilDiv X) where
  carrier := {D | IsCartierDiv X D}
  zero_mem' := isCartierDiv_zero X hnorm
  add_mem' hD hE := isCartierDiv_add X hD hE
  neg_mem' hD := isCartierDiv_neg X hD

/-! ## ★3. `Q`-Cartier 性を単系論の語彙へ(鎖 `cartier` の `qcartier`)

★`Found/Divisor/CartierMonoid.lean` の `IsQCartierSubgroup Γ := ∀ s, ∃ n>0, single s n ∈ Γ`
と、幾何側の「各素因子が `Q`-Cartier」が一致する、というのがこの節点である。 -/

set_option linter.unusedSectionVars false in
/-- ★★**各素因子が `Q`-Cartier ⟹ `Γ` は `Q`-Cartier 部分群**。

★`Γ` は `S`(素因子の部分集合)へ制限した Cartier 因子の群。
★これが `Found/Divisor/Cartier*.lean` の全定理の前提 `hQ` を与える。

★★★2026-09-06 に埋まった。**逸脱も追加仮定も無い** —— この 1 本だけは
`IsCartierDiv` の中身に一切触れず、`Finsupp.mapDomain_single` と `hΓ` の往復で閉じる。
★したがって `ordAtDiv` がまだ `sorry` だった時点で既に通っていた(本ファイルで唯一)。
残りの 3 本は同日 D8 の配線で埋まった。 -/
theorem isQCartierSubgroup_of_forall_isQCartier
    {S : Type u} (ι : S → PrimeDivisorPt X) (_hι : Function.Injective ι)
    (Γ : AddSubgroup (S →₀ ℤ))
    (_hΓ : ∀ x : S →₀ ℤ, x ∈ Γ ↔ IsCartierDiv X (Finsupp.mapDomain ι x))
    (_hQ : ∀ s : S, IsQCartierDiv X (Finsupp.single (ι s) 1)) :
    ABC3.Found.FrdI.IsQCartierSubgroup Γ := by
  intro s
  obtain ⟨n, hn, hc⟩ := _hQ s
  refine ⟨n, hn, (_hΓ (Finsupp.single s (n : ℤ))).mpr ?_⟩
  have hnat : (n • (Finsupp.single (ι s) (1 : ℤ)) : WeilDiv X)
      = Finsupp.single (ι s) (n : ℤ) := by
    rw [← Nat.cast_smul_eq_nsmul ℤ, Finsupp.smul_single, smul_eq_mul, mul_one]
  rw [Finsupp.mapDomain_single, ← hnat]
  exact hc

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def IsCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier divisor",
    sectionId := "frdi-example-6-1" }

def IsQCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Q-Cartier prime divisor",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_zero.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(零元)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_zero.needs : List ProofObligation :=
  [ .derivation "`f = 1` を取れば `div(1) = 0`" 109 ]

def isCartierDiv_add.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(和)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_add.needs : List ProofObligation :=
  [ .derivation "2 つの近傍の共通部分を取り、`f·g` を当てる(`ordAtDiv_mul`)" 109,
    .citation "[ABC3]" "ordAtDiv_mul" (.inProject "ABC3" "ABC3.Skeleton.Divisor.ordAtDiv_mul") 109 ]

def isCartierDiv_neg.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(逆元)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_neg.needs : List ProofObligation :=
  [ .derivation "`f⁻¹` を当てる" 109 ]

def cartierSubgroup.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Φ(L)^gp ⊆ ℤ[D_L] は Cartier 因子の群",
    sectionId := "frdi-example-6-1" }

def isQCartierSubgroup_of_forall_isQCartier.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — K-Q-Cartier ⟹ IsQCartierSubgroup",
    sectionId := "frdi-example-6-1" }

/-- ★★これが出れば `Found/Divisor/Cartier*.lean` の**全定理の前提が揃う**。 -/
def isQCartierSubgroup_of_forall_isQCartier.needs : List ProofObligation :=
  [ .citation "[ABC3]" "IsQCartierSubgroup(単系側の定義)"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsQCartierSubgroup") 109,
    .derivation "`n·[s]` が Cartier なら `Finsupp.single s n ∈ Γ`(`Finsupp.mapDomain_single` で `ι` を通す。★単射性 `hι` は使わなかった)" 109,
    .implicitStep "★原文は `K-Q-Cartier` を**仮定**として置く(「we shall assume that DK is K-Q-Cartier」)" 109 ]

end ABC3.Skeleton.Divisor
