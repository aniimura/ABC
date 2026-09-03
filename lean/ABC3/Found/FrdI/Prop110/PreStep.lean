/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import Mathlib.NumberTheory.PrimeCounting
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.MonoidPrime
import ABC3.Found.FrdI.Prop110.Vocabulary

/-!
# Prop110 —— (ii) の pre-step と Frobenius・`Div` の公式・(iii)

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★★(ii) —— **pre-step と Frobenius 型は順序を入れ替えられる**

原文 (FrdI p.34):
> (ii) Any composite morphism β ◦α of C, where α is a pre-step, and β is of

原文 (FrdI p.34):
> Frobenius type, may be written as a composite

★**3 本目の行(「where α′ is a pre-step, and β′ is of Frobenius type such that:」)は
引用できない**(事故 #3 の 5 度目。`′` が layout モードで消え 6/50 文字で停止)。
400 dpi 目視では読める。読みは下に日本語で書く。

★**合成の向き**: 原文の `β ◦ α` は「先に `α`、次に `β`」＝ Lean の `α ≫ β`。
`α′ ◦ β′` は「先に `β′`、次に `α′`」＝ `β′ ≫ α′`。
★**つまり「pre-step のあと Frobenius 型」を「Frobenius 型のあと pre-step」に組み替える。**

★原文の証明(p.35、目視)は
> assertion (ii) follows **formally** from assertion (i).
(直前に「in light of the existence of morphisms of Frobenius type of a given
Frobenius degree — cf. Definition 1.3, (ii)」と断ってある)。

★★**「formally」の中身は 3 手**:
1. `A` から次数 `degFr β` の Frobenius 型射 `β′` を取る(`frobDegSurj`)
2. **(i) の第2の場合**を `φ := α`(pre-step)と縦 `β′` に当てて、
   `α ≫ β'' = β′ ≫ α''`(`β''` は `X` から出る次数 `degFr β` の Frobenius 型)を得る
3. `β''` と `β` は **`X` から出る同次数の Frobenius 型射**なので `frobDegUniq` が同型 `ε` を与え、
   ★**`α′ := α'' ≫ ε`。**

★**原文が「formally」と言うのは正しい** —— 新しい材料は 1 つも要らず、
**(i) と `Definition 1.3, (ii)` だけで閉じる。**
★**ただし手数は 3 手ある。**「formally」は「短い」ではなく「新しい材料が要らない」の意味である。
-/

include P in
/-- **(ii)** —— `pre-step ≫ Frobenius 型` を `Frobenius 型 ≫ pre-step` に組み替える。

★次数は保たれる(`degFr β′ = degFr β`)。 -/
theorem prop_1_10_ii (F : FrobenioidCore P) {A X B : C}
    (α : A ⟶ X) (hα : IsPreStep P α) (β : X ⟶ B) (hβ : IsFrobeniusType P β) :
    ∃ (Y : C) (β' : A ⟶ Y) (α' : Y ⟶ B),
      IsFrobeniusType P β' ∧ P.degFr β' = P.degFr β ∧
      IsPreStep P α' ∧ β' ≫ α' = α ≫ β := by
  obtain ⟨Y, β', hβ', hdβ'⟩ := F.frobDegSurj A (P.degFr β)
  obtain ⟨X'', β'', α'', hβ'', hd'', hpre'', hsq⟩ :=
    prop_1_10_i_exists_preStep P F α hα β' hβ'
  -- `β''` と `β` は `X` から出る同次数の Frobenius 型射
  have hdeg : P.degFr β'' = P.degFr β := by rw [hd'', hdβ']
  obtain ⟨ε, hεiso, hε⟩ := F.frobDegUniq X X'' B β'' β hβ'' hβ hdeg
  haveI : IsIso ε := hεiso
  refine ⟨Y, β', α'' ≫ ε, hβ', hdβ',
    IsPreStep.comp P hpre'' (isPreStep_of_isIso P ε), ?_⟩
  rw [← Category.assoc, ← hsq, Category.assoc, hε]

/-! ### ★(ii) の `Div` の公式 —— **(i) のものがそのまま使える**

原文 (FrdI p.34):
> ∗for the bijection induced by applying the functor Φ to the base-

★**公式そのものの行(「Div(α′) = degFr(β) · β′∗(Div(α))」)は引用できない**
(事故 #3 の 6 度目、5/19 文字で停止)。隣接する `β′∗` の説明行を引く。

★★**(i) で作った `prop_1_10_i_Div_formula` が、四角形の向きを合わせるだけで
そのまま (ii) に効く。** 対応は
```
(i):  φ ≫ β_上 = α_左 ≫ φ′
(ii): α ≫ β    = β′   ≫ α′
```
で、`φ := α`(pre-step)、`β_上 := β`、`α_左 := β′`、`φ′ := α′`。

★**要る仮定は「`β′` と `β` が isometry」**で、どちらも Frobenius 型だから出る。
★**(i) で「原文が書いていない材料」として明示した isometry が、
(ii) でも同じ位置に要る** —— つまりあれは (i) 固有の抜けではなく、
**この形の公式に共通して要るもの**だった。

★**実体は下(`prop_1_10_ii_Div_formula`)に置く** —— (i) の公式を使うので、
**宣言の順序でそちらが先に要る**。★これは Lean の制約であって数学の順序ではない。
-/

include P in
/-- **(i) の `Div` の公式**(`Φ.map (Base α)` を当てた形)。

★**`α` と `β` が isometry であることを仮定に明示している** ——
原文は「Remark 1.1.1 から」と書くが、**実際にはこれが要る**。
★仮定を文に出したので、**何が効いているかが宣言の型から読める。** -/
theorem prop_1_10_i_Div_formula {A B A' B' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (β : B ⟶ B') (φ' : A' ⟶ B')
    (hα : IsIsometric P α) (hβ : IsIsometric P β)
    (hsq : φ ≫ β = α ≫ φ') :
    Φ.map (P.Base α) (P.Div φ') = (P.degFr β : ℕ) • P.Div φ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hα, hβ] at h
  simpa using h.symm

include P in
/-- ★**(i) の `Div` の公式 —— 原文の `α∗`(全単射)の向き**。

★上の `prop_1_10_i_Div_formula` は `Φ.map (Base α)` を**左辺に押し出した**形で、
原文の `Div(φ′) = degFr(β) · α∗(Div(φ))` の形ではない。
`α` は Frobenius 型ゆえ base-isomorphism なので、その逆で向きを直せる。

★★**検証役の指摘**: (ii) の同形の補題(`prop_1_10_ii_Div_formula'`)が
**(ii) の名前の下にしかなく、(i) を探した読者は押し出し形しか見つけられない**。
★**変数の対応は (ii) の `(α, β, β′, α′)` に (i) の `(φ, β, α, φ′)` を代入したものだが、
名前が違えば見つからない。** -/
theorem prop_1_10_i_Div_formula' {A B A' B' : C}
    (φ : A ⟶ B) (α : A ⟶ A') (β : B ⟶ B') (φ' : A' ⟶ B')
    (hα : IsFrobeniusType P α) (hβ : IsIsometric P β)
    (hsq : φ ≫ β = α ≫ φ') :
    haveI : IsIso (P.Base α) := hα.2
    P.Div φ' = (P.degFr β : ℕ) • Φ.map (inv (P.Base α)) (P.Div φ) := by
  haveI : IsIso (P.Base α) := hα.2
  have h := prop_1_10_i_Div_formula P φ α β φ' hα.1.2 hβ hsq
  have h2 := congrArg (Φ.map (inv (P.Base α))) h
  rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at h2
  exact h2

/-- **(i) の `degFr` の等式**(第1の場合の構成に対して)。

★上の `prop_1_10_i_exists_frobType` が返す `φ′` は `ψ ≫ γ` の形で、
`degFr ψ = degFr φ`、`γ` は同型。★**だから次数は保たれる。**

★**`FrobenioidCore` は要らない** —— `degFr` の関手性と同型の次数だけで出る。
★**「原文が (i) の材料として `Definition 1.3, (ii)` を挙げている」のは
存在のためであって、次数の等式のためではない。** -/
theorem prop_1_10_i_degFr_eq {A' Z B' : C}
    (ψ : A' ⟶ Z) (γ : Z ⟶ B') [IsIso γ] {n : ℕ+} (hψ : P.degFr ψ = n) :
    P.degFr (ψ ≫ γ) = n := by
  rw [P.degFr_comp, degFr_of_isIso, one_mul, hψ]

include P in
/-- **(ii) の `Div` の公式**。★(i) のものを向きを合わせて使うだけ。

★★**(i) で「原文が書いていない材料」として明示した isometry が、
(ii) でも同じ位置に要る** —— つまりあれは (i) 固有の抜けではなく、
**この形の公式に共通して要るもの**だった。 -/
theorem prop_1_10_ii_Div_formula {A X B Y : C}
    (α : A ⟶ X) (β : X ⟶ B) (β' : A ⟶ Y) (α' : Y ⟶ B)
    (hβ' : IsFrobeniusType P β') (hβ : IsFrobeniusType P β)
    (hsq : β' ≫ α' = α ≫ β) :
    Φ.map (P.Base β') (P.Div α') = (P.degFr β : ℕ) • P.Div α :=
  prop_1_10_i_Div_formula P α β' β α' hβ'.1.2 hβ.1.2 hsq.symm

include P in
/-- ★★★**原文の形の `Div` の公式**（2026-08-16）。

原文 (FrdI p.34):
> ∗for the bijection induced by applying the functor Φ to the base-

★**引用を選び直した記録（事故 #3 の 9 度目）**: 主張の行
「where α′ is a pre-step, and β′ is of Frobenius type such that:」は
★**`′` が抽出で落ちる**ため引用できない（6/50 文字で停止）。
★`′` は本ファイルで何度も落ちている既知の文字である。

★★**監査で「`β′` の base-isomorphism 性を仮定せず、原文の `β′∗`（全単射）の
形になっていない」と指摘されたもの**である。

★`β′` は Frobenius 型なので base-isomorphism、したがって
`Φ.map (Base β′)` は全単射であり、その逆を使って原文の向きに直せる。 -/
theorem prop_1_10_ii_Div_formula' {A X B Y : C}
    (α : A ⟶ X) (β : X ⟶ B) (β' : A ⟶ Y) (α' : Y ⟶ B)
    (hβ' : IsFrobeniusType P β') (hβ : IsFrobeniusType P β)
    (hsq : β' ≫ α' = α ≫ β) :
    haveI : IsIso (P.Base β') := hβ'.2
    P.Div α' = (P.degFr β : ℕ) • Φ.map (inv (P.Base β')) (P.Div α) := by
  haveI : IsIso (P.Base β') := hβ'.2
  have h := prop_1_10_ii_Div_formula P α β β' α' hβ' hβ hsq
  have h2 := congrArg (Φ.map (inv (P.Base β'))) h
  rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id, map_nsmul] at h2
  exact h2


/-! ## ★(iii) —— `of perfect type` と像のモノイドの perfect 性

原文 (FrdI p.34):
> (iii) Suppose that C is of perfect type. Then the monoids in the image of

原文 (FrdI p.34):
> Φ are perfect. If, moreover, C is of isotropic and Frobenius-normalized type,

★**原文の証明**(p.35、目視):
> In light of the existence of morphisms of Frobenius type of a given Frobenius degree
> [cf. Definition 1.3, (ii)] and the equivalences of categories of Definition 1.3, (iii), (d),
> the fact that Φ(A) is perfect follows immediately [cf. Remark 1.1.1] from the fact
> that A is perfect [cf. Definition 1.2, (iv)].

★**不足していた語彙 `of perfect type` をここで足す** ——
`IsPerfectObj` は対象について定義済みで、「型」は**全対象がそれである**という形。
-/

/-- **[FrdI] Definition 1.3 系** `𝒞` が **of perfect type** —— 全対象が perfect。

★**「型」の定義はどれもこの形**(全対象がその性質を持つ)なので、
`IsPerfectObj` から機械的に作れる。★**原文が別に定義を置いていないのは
そのためだが、我々は名前が要る。** -/
def IsOfPerfectType : Prop := ∀ A : C, IsPerfectObj P A

/-! ### ★★(iii) の核心 —— `n` 倍が全射になる仕組み

★**原文は「follows immediately [cf. Remark 1.1.1]」と書く。**
その `Remark 1.1.1`(＝ `Div_comp`)を四角形
`φ₁ ≫ ψ' = ψ ≫ φ₂`(`φ₁`・`φ₂` は次数 `n` の Frobenius 型、`ψ`・`ψ'` は pre-step)
に当てると:
```
Div(φ₁ ≫ ψ') = Φ.map (Base φ₁) (Div ψ') + (degFr ψ') • Div φ₁
Div(ψ ≫ φ₂)  = Φ.map (Base ψ)  (Div φ₂) + (degFr φ₂)  • Div ψ
```
★**`Div φ₁ = Div φ₂ = 0`**(Frobenius 型 ⟹ isometry)なので、両辺は
```
Φ.map (Base φ₁) (Div ψ') = n • Div ψ
```
★★**これが「`n` 倍が全射」の正体である** ——
`IsPerfectObj` の条件 (b) が `ψ` の存在を与え、その `Div ψ` が
`Φ.map (Base φ₁) (Div ψ')` の `n` 分の 1 になる。

★**(i) の `Div` の公式とまったく同じ形**である(`prop_1_10_i_Div_formula`)。
★**同じ補題が (i)・(ii)・(iii) の 3 箇所で効いている** ——
原文が 3 箇所で別々に「Remark 1.1.1 から」と書いているものは、**1 本の等式**だった。
-/

include P in
/-- **(iii) の核心** —— 次数 `n` の Frobenius 型射に沿った四角形は、
`Div` を **`n` 倍の関係**に翻訳する。

★`prop_1_10_i_Div_formula` の直接の帰結。★**新しい証明は 1 行も要らない。** -/
theorem prop_1_10_iii_div_nsmul {B₁ B₁' B₂ B₂' : C}
    (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂') (ψ : B₁ ⟶ B₂) (ψ' : B₁' ⟶ B₂')
    (hφ₁ : IsFrobeniusType P φ₁) (hφ₂ : IsFrobeniusType P φ₂)
    {n : ℕ+} (hn : P.degFr φ₂ = n)
    (hsq : φ₁ ≫ ψ' = ψ ≫ φ₂) :
    Φ.map (P.Base ψ) (P.Div φ₂) + (n : ℕ) • P.Div ψ
      = Φ.map (P.Base φ₁) (P.Div ψ') + (P.degFr ψ' : ℕ) • P.Div φ₁ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hn] at h
  exact h.symm

include P in
/-- ★**上を、両辺の消える項を落とした形で**。★**これが原文の言う
「Remark 1.1.1 から従う」の実体**である。 -/
theorem prop_1_10_iii_div_nsmul' {B₁ B₁' B₂ B₂' : C}
    (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂') (ψ : B₁ ⟶ B₂) (ψ' : B₁' ⟶ B₂')
    (hφ₁ : IsIsometric P φ₁) (hφ₂ : IsIsometric P φ₂)
    {n : ℕ+} (hn : P.degFr φ₂ = n)
    (hsq : φ₁ ≫ ψ' = ψ ≫ φ₂) :
    Φ.map (P.Base φ₁) (P.Div ψ') = (n : ℕ) • P.Div ψ := by
  have h := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, hn, hφ₁, hφ₂] at h
  simpa using h

/-! ## ★(v) —— Frobenius 型 ⟺ prime-Frobenius の合成

原文 (FrdI p.35):
> (v) A morphism of C is a morphism of Frobenius type if and only if it is a

原文 (FrdI p.35):
> composite of prime-Frobenius morphisms.

★★**「合成」を帰納的に定義するとき、原文の曖昧さが 2 段階で出る。**

**第1段**: 原文は「a composite of prime-Frobenius morphisms」と書くだけで、
**空の合成を含むかを言っていない**。含めないと、次数 1 の Frobenius 型射が
合成として書けない。

**第2段(こちらが本質)**: ★**空の合成を「恒等射」に取るだけでは足りない。**
`Proposition 1.4, (iii)` により **次数 1 の Frobenius 型射は同型**であり、
同型は一般に **`𝟙 A` そのものではない**。(v) が iff である以上、
★**基底の場合は「同型」でなければならない。**

★★**これが原文の 3 つ目の曖昧さである。**
「composite of prime-Frobenius morphisms」は **同型を法として読む**しかない。
★**形式化しなければ「恒等射」と「同型」の差は見えない** ——
自然言語では「合成が無い場合」で済んでしまう。
-/

/-- **prime-Frobenius 射の合成**(帰納的)。

★**基底は同型**(恒等射ではない)。理由は上の docstring を見よ ——
次数 1 の Frobenius 型射は同型であって、`𝟙` とは限らない。 -/
inductive IsPrimeFrobComposite : ∀ {A B : C}, (A ⟶ B) → Prop
  | iso {A B : C} (φ : A ⟶ B) : IsIso φ → IsPrimeFrobComposite φ
  | cons {A B E : C} {φ : A ⟶ B} {ψ : B ⟶ E} :
      IsPrimeFrobenius P φ → IsPrimeFrobComposite ψ → IsPrimeFrobComposite (φ ≫ ψ)

include P in
/-- **(v) の片向き** —— prime-Frobenius の合成は Frobenius 型。

★**`Proposition 1.7, (i)` の閉性だけで出る。** 原文が (v) の証明で
何を使うかは書いていないが、★**この向きは合成閉性の帰納法である。** -/
theorem isFrobeniusType_of_isPrimeFrobComposite (F : FrobenioidCore P)
    {A B : C} {φ : A ⟶ B} (h : IsPrimeFrobComposite P φ) : IsFrobeniusType P φ := by
  induction h with
  | iso ψ hψ => exact @isFrobeniusType_of_isIso _ _ _ _ _ P _ _ ψ hψ
  | cons hφ _ ih => exact IsFrobeniusType.comp P F hφ.1 ih

include P in
/-- ★**非退化**: prime-Frobenius 自身は(恒等射との合成として)合成である。
★**「合成の集合が空でない」を押さえる。** -/
theorem isPrimeFrobComposite_of_isPrimeFrobenius {A B : C} {φ : A ⟶ B}
    (h : IsPrimeFrobenius P φ) : IsPrimeFrobComposite P φ := by
  have := IsPrimeFrobComposite.cons h
    (IsPrimeFrobComposite.iso (P := P) (𝟙 B) inferInstance)
  simpa using this

/-! ### ★(v) の逆向き —— Frobenius 型 ⟹ prime-Frobenius の合成

★**原文は証明を書いていない**((v) の証明は本文に無い)。
★**構成は次数についての強帰納法**である:
- `degFr φ = 1` のとき: ★**`Proposition 1.4, (iii)` により `φ` は同型**。基底の場合。
- `degFr φ = n > 1` のとき: `n` の素因数 `p` を取り `n = m * p` と書く。
  次数 `p` の Frobenius 型射 `φ₁` と次数 `m` の `φ₂` を `frobDegSurj` で取ると
  `φ₁ ≫ φ₂` の次数は `m * p = n` なので、`frobDegUniq` が同型 `δ` を与え
  ★**`φ = φ₁ ≫ (φ₂ ≫ δ)`**。`φ₂ ≫ δ` の次数は `m < n` なので帰納法が回る。

★★**「素因数を 1 つ剥がす」が帰納のステップである。** 原文の (v) は
「Frobenius 型 ⟺ prime-Frobenius の合成」と書くだけだが、
★**その中身は「次数の素因数分解」そのもの**である。
-/

include P in
/-- **(v) の逆向き** —— Frobenius 型射は prime-Frobenius 射の合成。

★次数についての強帰納法。基底は **同型**(次数 1)。 -/
theorem isPrimeFrobComposite_of_isFrobeniusType (F : FrobenioidCore P) :
    ∀ (n : ℕ) {A B : C} (φ : A ⟶ B), IsFrobeniusType P φ →
      ((P.degFr φ : ℕ+) : ℕ) = n → IsPrimeFrobComposite P φ := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro A B φ hφ hn
    rcases eq_or_ne n 1 with h1 | h1
    · -- 次数 1 ⟹ 同型
      have hlin : IsLinear P φ := by
        show P.degFr φ = 1
        exact PNat.coe_injective (by rw [hn, h1]; rfl)
      exact IsPrimeFrobComposite.iso φ (prop_1_4_iii P F φ hφ.1 ⟨hlin, hφ.2⟩)
    · -- 素因数を 1 つ剥がす
      obtain ⟨p, hp, k, hk⟩ := Nat.exists_prime_and_dvd h1
      have hnpos : 0 < n := by rw [← hn]; exact (P.degFr φ).pos
      have hppos : 0 < p := hp.pos
      have hkpos : 0 < k := Nat.pos_of_ne_zero (by rintro rfl; omega)
      set pp : ℕ+ := ⟨p, hppos⟩ with hpp
      set kk : ℕ+ := ⟨k, hkpos⟩ with hkk
      obtain ⟨Z, φ₁, hφ₁, hd₁⟩ := F.frobDegSurj A pp
      obtain ⟨W, φ₂, hφ₂, hd₂⟩ := F.frobDegSurj Z kk
      have hcomp : IsFrobeniusType P (φ₁ ≫ φ₂) := IsFrobeniusType.comp P F hφ₁ hφ₂
      have hdc : P.degFr (φ₁ ≫ φ₂) = P.degFr φ := by
        rw [P.degFr_comp, hd₁, hd₂]
        refine PNat.coe_injective ?_
        simp only [PNat.mul_coe, hpp, hkk, hn, hk]
        exact Nat.mul_comm k p
      obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A W B (φ₁ ≫ φ₂) φ hcomp hφ hdc
      haveI : IsIso δ := hδiso
      have hdk : ((P.degFr (φ₂ ≫ δ) : ℕ+) : ℕ) = k := by
        rw [P.degFr_comp, degFr_of_isIso P δ, one_mul, hd₂, hkk]
        rfl
      have hklt : k < n := by
        rw [hk]
        calc k = 1 * k := (Nat.one_mul k).symm
          _ < p * k := Nat.mul_lt_mul_of_lt_of_le hp.one_lt (le_refl k) hkpos
      have hrec := ih k hklt (φ₂ ≫ δ)
        (IsFrobeniusType.comp P F hφ₂ (isFrobeniusType_of_isIso P δ)) hdk
      have hprime : IsPrimeFrobenius P φ₁ := ⟨hφ₁, by rw [hd₁]; simpa [hpp] using hp⟩
      have : φ = φ₁ ≫ φ₂ ≫ δ := by rw [← hδ, Category.assoc]
      rw [this]
      exact IsPrimeFrobComposite.cons hprime hrec

include P in
/-- ★★**`Proposition 1.10, (v)` の完成形** —— Frobenius 型 ⟺ prime-Frobenius の合成。

★**原文はこの iff を証明なしで述べる。** 我々の側では
- `⟸` は `Proposition 1.7, (i)` の合成閉性の帰納法
- `⟹` は **次数の素因数分解による強帰納法**(基底は `Proposition 1.4, (iii)` で同型)
であり、★**「合成」の基底を同型に取らないと iff にならない**ことが分かった。 -/
theorem prop_1_10_v (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsFrobeniusType P φ ↔ IsPrimeFrobComposite P φ :=
  ⟨fun h => isPrimeFrobComposite_of_isFrobeniusType P F _ φ h rfl,
   isFrobeniusType_of_isPrimeFrobComposite P F⟩

/-! ### ★★(v) の基底 —— 「意図的な修復」を**反例**に格上げする(2026-08-16)

★★**取り下げ表の (v) の項目**: 基底を「任意の同型」に取ったことを
**意図的な修復**として記録していた。★しかし「修復」と「反例がある」は違う ——
★**素直な読み(基底は恒等射)が実際に偽であることを示せば、
同型に取るのは選択ではなく強制になる。**

★**ここでそれを示す。** 素直な読み `IsPrimeFrobCompositeId` を別に定義し、
★**恒等射でない同型は決してその形に書けない**ことを証明する。
`Proposition 1.7, (i)`(同型は Frobenius 型)により、そのような射は
**(v) の左辺を満たす**から、素直な読みでは iff が破れる。

★★**理由は次数だけである** —— prime-Frobenius 射の次数は素数(≥2)であり、
次数は合成で掛かる(`Remark 1.1.1`)。★空でない合成の次数は 2 以上だが、
同型の次数は 1 である。★**だから空の合成しか残らず、それは `𝟙` である。**

★★**これは「原文の言葉が形式化で 2 通りに分かれ、一方が偽になる」例である。**
自然言語の「composite of prime-Frobenius morphisms」は、
★**同型を法として読む以外にない。**
-/

/-- ★**(v) の素直な読み** —— 基底を**恒等射**に取った「prime-Frobenius の合成」。

★これが原文の字面どおりの読みである。★下でこれが**偽**であることを示す。 -/
inductive IsPrimeFrobCompositeId : ∀ {A B : C}, (A ⟶ B) → Prop
  | nil {A : C} : IsPrimeFrobCompositeId (𝟙 A)
  | cons {A B E : C} {φ : A ⟶ B} {ψ : B ⟶ E} :
      IsPrimeFrobenius P φ → IsPrimeFrobCompositeId ψ → IsPrimeFrobCompositeId (φ ≫ ψ)

include P in
/-- ★**次数 1 の「素直な合成」は空の合成しかない**。

★prime-Frobenius 射の次数は素数(≥2)で、次数は合成で掛かるので、
★**`cons` の場合は次数が 2 以上になり、1 にはなり得ない。** -/
theorem eq_id_of_isPrimeFrobCompositeId_of_degFr_one :
    ∀ {A B : C} {φ : A ⟶ B}, IsPrimeFrobCompositeId P φ →
      ((P.degFr φ : ℕ+) : ℕ) = 1 → ∃ h : A = B, φ = eqToHom h := by
  intro A B φ h
  induction h with
  | nil => exact fun _ => ⟨rfl, by simp⟩
  | @cons A' B' E' φ' ψ' hp _ _ =>
      intro hd
      exfalso
      rw [P.degFr_comp, PNat.mul_coe] at hd
      have h2 : 2 ≤ ((P.degFr φ' : ℕ+) : ℕ) := hp.2.two_le
      have h1 : ((P.degFr φ' : ℕ+) : ℕ) = 1 :=
        Nat.dvd_one.mp (Dvd.intro_left _ hd)
      omega

include P in
/-- ★★★**(v) の素直な読みは偽である** —— 恒等射でない自己同型は、
prime-Frobenius 射の(恒等射を基底とする)合成には決して書けない。

★★**しかし `Proposition 1.7, (i)` によりそれは Frobenius 型である。**
したがって素直な読みでは (v) の iff が破れる ——
★**基底を同型に取ったのは修復ではなく、強制であった。** -/
theorem not_isPrimeFrobCompositeId_of_isIso_of_ne_id {A : C} (u : A ⟶ A) [IsIso u]
    (hne : u ≠ 𝟙 A) : ¬ IsPrimeFrobCompositeId P u := by
  intro h
  obtain ⟨hAA, hu⟩ := eq_id_of_isPrimeFrobCompositeId_of_degFr_one P h
    (by rw [degFr_of_isIso P u]; rfl)
  exact hne (by simpa using hu)

/-! ## ★(vi) —— `𝒞^istr` と group-like 対象

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

原文 (FrdI p.35):
> every group-like object A ∈Ob(Cistr) is Frobenius-trivial.

★**後半(group-like ⟹ Frobenius-trivial)の証明**(p.36、目視):
> If, moreover, A is group-like, then [since Cistr is a Frobenioid — cf. Proposition 1.9, (v)]
> it follows from Definition 1.3, (i), (a), (b), that there exist [co-angular — cf.
> Proposition 1.4, (i)] pre-steps A′ →A, A′ →A′′, where A′′ is Frobenius-trivial.
> But by Proposition 1.4, (iii), these pre-steps are isomorphisms, so A is Frobenius-trivial.

★★**「これらの pre-step が同型である」に `Proposition 1.4, (iii)` を使うには
LB-invertible が要る。** 原文は `[co-angular — cf. Proposition 1.4, (i)]` で
**co-angular のほうだけ**括弧に入れており、★**isometric のほうを書いていない。**

★**isometric はどこから来るか**: `A` が group-like で、`Φ` が **divisorial(＝ sharp を含む)**
なら、★**`Φ(A)` のすべての元が 0** である。よって `A` から出る任意の射の `Div` は 0、
すなわち **isometric**。

★★**これは (i) の `Div` の公式・(ii)・(iii) で見つけた「原文が isometry を書かない」
というパターンの 4 例目**である。★**同じ省略が 4 箇所で起きている。**
-/

/-- ★**group-like ＋ sharp ⟹ モノイドは零**。

★`IsGroupLike M ↔ ∀ a, IsAddUnit a`(`isGroupLike_iff`)と
`IsSharp M := ∀ a, IsAddUnit a → a = 0` を合わせるだけ。
★**しかしこの一行が (vi) の後半の要点である。** -/
theorem eq_zero_of_isGroupLike_of_isSharp {M : Type w} [AddCommMonoid M]
    (hg : IsGroupLike M) (hs : IsSharp M) (a : M) : a = 0 :=
  hs a ((isGroupLike_iff M).mp hg a)

include P in
/-- ★★**group-like 対象から出る射はすべて isometric**。

★**原文が書いていない一歩**(上の docstring を見よ)。
`Φ` が divisorial であることの **sharp の部分**がここで効く ——
★**divisorial の定義に sharp が入っている理由の 1 つがこれである。** -/
theorem isIsometric_of_isGroupLikeObj {A B : C} (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) : IsIsometric P φ :=
  eq_zero_of_isGroupLike_of_isSharp hA (P.divisorial _).2 (P.Div φ)

include P in
/-- ★**(vi) の後半の要点** —— isotropic かつ group-like な対象から出る
pre-step は同型。

★原文の「by Proposition 1.4, (iii), these pre-steps are isomorphisms」の中身。
★**co-angular は isotropic から(`Proposition 1.4, (i)`)、
isometric は group-like から(上の補題)。原文は後者を書いていない。** -/
theorem isIso_of_preStep_of_isGroupLikeObj (F : FrobenioidCore P) {A B : C}
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X) (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) (hφ : IsPreStep P φ) : IsIso φ :=
  prop_1_4_iii P F φ ⟨prop_1_4_i P φ hiso, isIsometric_of_isGroupLikeObj P hA φ⟩ hφ

include P in
/-- ★**`𝒞` が isotropic 型のときの形**(`𝒞^istr` で使うのはこちら)。

★`Proposition 1.9` の `𝒞^istr` は全対象が isotropic なので、仮定はそこから供給される。 -/
theorem isIso_of_preStep_of_isGroupLikeObj' (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (hA : IsGroupLikeObj P A)
    (φ : A ⟶ B) (hφ : IsPreStep P φ) : IsIso φ :=
  isIso_of_preStep_of_isGroupLikeObj P F (fun X _ => hiso X) hA φ hφ

/-! ## ★(iv) —— prime-Frobenius ⟺ irreducible

原文 (FrdI p.34):
> nius morphism if and only if it is irreducible [cf. §0]. In particular, if A ∈Ob(C)

★**原文の証明**(p.36、目視):
> assertion (iv) follows immediately from Proposition 1.7, (v), and
> **the well-known structure of the multiplicative monoid N≥1**
> [cf. also Definition 1.3, (ii)]

★★**「the well-known structure of the multiplicative monoid ℕ≥1」とは何か。**
原文はこの一句で (iv) を片づける。★**その中身は「乗法モノイド `(ℕ≥1, ×)` において
既約な元とは素数のことである」**である。

★**mathlib には `PNat.Prime`(＝ `(p : ℕ).Prime`)はあるが、
「乗法モノイド `ℕ+` の既約元」という形の補題は無い**(`Mathlib/Data/PNat/Prime.lean` を確認)。
★**だからここで書く。** ——★**原文が「well-known」と呼ぶものが、
我々の道具立てには存在しなかった、という測定でもある。**

★**そしてこれが (iv) の骨格である**: 射 `φ` の既約性は
「`φ = β ≫ α` の分解で片方が同型」であり、Frobenius 型射に限れば
**次数の分解 `degFr φ = degFr α * degFr β`** に写る。同型は次数 1 なので、
★**射の既約性は次数の既約性に、そのまま翻訳される。**
-/

/-- ★★**乗法モノイド `ℕ≥1` の既約元は素数**(原文の「well-known structure」)。

★`n : ℕ+` について「1 でなく、どう 2 つに分けても片方が 1」⟺「素数」。
★**mathlib に無いので書いた。** -/
theorem pnat_irreducible_iff_prime (n : ℕ+) :
    (n ≠ 1 ∧ ∀ a b : ℕ+, n = a * b → a = 1 ∨ b = 1) ↔ Nat.Prime (n : ℕ) := by
  have hn0 : (n : ℕ) ≠ 0 := n.pos.ne'
  constructor
  · rintro ⟨hne, hsplit⟩
    rw [← Nat.irreducible_iff_nat_prime]
    refine ⟨fun h => hne (PNat.coe_injective ?_), fun a b hab => ?_⟩
    · rw [Nat.isUnit_iff] at h; rw [h]; rfl
    · have ha0 : a ≠ 0 := by rintro rfl; rw [Nat.zero_mul] at hab; exact hn0 hab
      have hb0 : b ≠ 0 := by rintro rfl; rw [Nat.mul_zero] at hab; exact hn0 hab
      have hsp := hsplit ⟨a, Nat.pos_of_ne_zero ha0⟩ ⟨b, Nat.pos_of_ne_zero hb0⟩
        (PNat.coe_injective (by simpa using hab))
      rcases hsp with h | h
      · exact Or.inl (Nat.isUnit_iff.mpr (congrArg PNat.val h))
      · exact Or.inr (Nat.isUnit_iff.mpr (congrArg PNat.val h))
  · intro hp
    have hirr : Irreducible (n : ℕ) := (Nat.irreducible_iff_nat_prime _).mpr hp
    refine ⟨fun h => ?_, fun a b hab => ?_⟩
    · rw [h] at hp; exact Nat.not_prime_one (by simpa using hp)
    · have hab' : (n : ℕ) = (a : ℕ) * (b : ℕ) := by rw [hab]; rfl
      rcases hirr.isUnit_or_isUnit hab' with h | h
      · exact Or.inl (PNat.coe_injective (by simpa using Nat.isUnit_iff.mp h))
      · exact Or.inr (PNat.coe_injective (by simpa using Nat.isUnit_iff.mp h))

/-! ### ★(iv) の片向き —— irreducible ⟹ prime-Frobenius

★**射の分解が次数の分解に写る**のが要点:
`φ = φ₁ ≫ (φ₂ ≫ δ)` の形を `frobDegSurj` ＋ `frobDegUniq` で作り、
既約性から片方が同型と分かると、★**その次数が 1 になる**。
つまり **`degFr φ = a * b` ⟹ `a = 1` または `b = 1`**。
そこに `pnat_irreducible_iff_prime` を当てて素数性を得る。

★**`isotropic` の仮定が効く場所**: 「次数 1 ⟹ 同型」に `Proposition 1.4, (iii)` を
使うには LB-invertible が要り、その co-angular の部分を
**`Proposition 1.4, (i)`(isotropic なら全射が co-angular)**が供給する。
★**原文が `[cf. §0]` としか書かない (iv) の中で、isotropic はここに効いている。**
-/

include P in
/-- **(iv) の片向き** —— isotropic な域を持つ Frobenius 型射が irreducible なら
prime-Frobenius。 -/
theorem prop_1_10_iv_mpr (F : FrobenioidCore P) {A B : C}
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X)
    (φ : A ⟶ B) (hφ : IsFrobeniusType P φ) (hirr : IsIrreducibleMor φ) :
    IsPrimeFrobenius P φ := by
  refine ⟨hφ, ?_⟩
  rw [← pnat_irreducible_iff_prime]
  refine ⟨fun h => hirr.1 (prop_1_4_iii P F φ hφ.1 ⟨h, hφ.2⟩), fun a b hab => ?_⟩
  obtain ⟨Z, φ₁, hφ₁, hd₁⟩ := F.frobDegSurj A a
  obtain ⟨W, φ₂, hφ₂, hd₂⟩ := F.frobDegSurj Z b
  have hcomp : IsFrobeniusType P (φ₁ ≫ φ₂) := IsFrobeniusType.comp P F hφ₁ hφ₂
  have hdc : P.degFr (φ₁ ≫ φ₂) = P.degFr φ := by
    rw [P.degFr_comp, hd₁, hd₂, hab, mul_comm]
  obtain ⟨δ, hδiso, hδ⟩ := F.frobDegUniq A W B (φ₁ ≫ φ₂) φ hcomp hφ hdc
  haveI : IsIso δ := hδiso
  have hfac : φ = φ₁ ≫ (φ₂ ≫ δ) := by rw [← hδ, Category.assoc]
  rcases hirr.2 Z φ₁ (φ₂ ≫ δ) hfac with h | h
  · haveI := h
    exact Or.inl (by rw [← hd₁]; exact (degFr_of_isIso P φ₁).symm ▸ rfl)
  · haveI := h
    refine Or.inr ?_
    have hd := degFr_of_isIso P (φ₂ ≫ δ)
    rw [P.degFr_comp, degFr_of_isIso P δ, one_mul, hd₂] at hd
    exact hd

end ABC3.Found.FrdI
