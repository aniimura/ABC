/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import Mathlib.NumberTheory.PrimeCounting
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.MonoidPrime
import ABC3.Found.FrdI.Prop110.Moreover

/-!
# Prop110 —— (iv) の In particular・(vi) 後半・Frobenius-compact

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)
variable {M N : Type*} [Monoid M] [Monoid N]

/-! ### ★(iv) の「In particular」—— 無限個の同型類

原文 (FrdI p.34):
> that arise from irreducible arrows with domain A.

★**構成**: 各素数 `p` に対し `frobDegSurj` が次数 `p` の Frobenius 型射 `φ_p : A ⟶ B_p` を与え、
(iv) によりそれは **irreducible**。
★**異なる素数は非同型な対象を与える**: `_A𝒞` の同型 `θ : B_p ⟶ B_q`(`φ_p ≫ θ = φ_q`)が
あれば `degFr φ_q = degFr θ * degFr φ_p = 1 * p = p`。★**`p ≠ q` に矛盾。**

★**素数が無限にあることは `Nat.exists_infinite_primes`(mathlib)から。**
★**「無限個の同型類」は「素数から同型類への単射がある」として書く** ——
原文の "infinitely many isomorphism classes" の忠実な形である。
-/

include P in
/-- ★**異なる次数の Frobenius 型射は `_A𝒞` で非同型**。

★これが「無限個の同型類」の核心。次数が同型類の不変量になっている。 -/
theorem frobType_not_iso_of_degFr_ne {A Bp Bq : C}
    (φp : A ⟶ Bp) (φq : A ⟶ Bq) (hne : P.degFr φp ≠ P.degFr φq)
    (θ : Bp ⟶ Bq) [IsIso θ] : φp ≫ θ ≠ φq := by
  intro h
  apply hne
  have : P.degFr (φp ≫ θ) = P.degFr φq := by rw [h]
  rwa [P.degFr_comp, degFr_of_isIso P θ, one_mul] at this

include P in
/-- ★★**(iv) の「In particular」** —— isotropic な `A` について、
**各素数ごとに irreducible な射があり、異なる素数のものは非同型**。

★**これが「無限個の同型類」の内容である**(素数は無限にあるので)。 -/
theorem prop_1_10_iv_infinitely_many (F : FrobenioidCore P)
    {A : C} (hA : IsIsotropic P A) :
    (∀ p : ℕ+, Nat.Prime (p : ℕ) →
       ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ P.degFr φ = p) ∧
    (∀ (Bp Bq : C) (φp : A ⟶ Bp) (φq : A ⟶ Bq), P.degFr φp ≠ P.degFr φq →
       ∀ θ : Bp ≅ Bq, φp ≫ θ.hom ≠ φq) := by
  constructor
  · intro p hp
    obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A p
    refine ⟨B, φ, ?_, hd⟩
    exact prop_1_10_iv_mp P F hA φ ⟨hφ, by rw [hd]; exact hp⟩
  · intro Bp Bq φp φq hne θ
    haveI : IsIso θ.hom := inferInstance
    exact frobType_not_iso_of_degFr_ne P φp φq hne θ.hom

include P in
/-- ★★**原文「there exist infinitely many isomorphism classes」そのもの**。

原文 (FrdI p.34):
> is isotropic, then there exist infinitely many isomorphism classes of objects of AC

★**監査以前は「無限個」を述べていなかった**(素数ごとの存在と非同型性まで)。
★**次数がいくらでも大きく取れる**という形で述べる ——
`prop_1_10_iv_infinitely_many` の第2主張(次数が違えば非同型)と合わせて、
★同型類が無限個あることを与える。 -/
theorem prop_1_10_iv_unbounded (F : FrobenioidCore P) {A : C} (hA : IsIsotropic P A) :
    ∀ n : ℕ, ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ n ≤ (P.degFr φ : ℕ) := by
  intro n
  obtain ⟨p, hpn, hp⟩ := Nat.exists_infinite_primes n
  have hppos : 0 < p := hp.pos
  obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A ⟨p, hppos⟩
  refine ⟨B, φ, prop_1_10_iv_mp P F hA φ ⟨hφ, ?_⟩, ?_⟩
  · rw [hd]; exact hp
  · rw [hd]; exact hpn

include P in
/-- ★★★**原文の「infinitely many isomorphism classes」を集合の無限性で述べる**。

原文 (FrdI p.34):
> is isotropic, then there exist infinitely many isomorphism classes of objects of AC

★★**監査で「「無限個」そのものを述べていない」と指摘されたもの**。

★**実現される次数の集合が無限**であることを述べる。
★`prop_1_10_iv_infinitely_many` の第2主張（次数が違えば非同型）と合わせれば、
★★**同型類が無限個ある**ことになる（次数が全単射の代わりをする）。 -/
theorem prop_1_10_iv_degrees_infinite (F : FrobenioidCore P) {A : C}
    (hA : IsIsotropic P A) :
    {n : ℕ | ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧ ((P.degFr φ : ℕ+) : ℕ) = n}.Infinite := by
  refine Set.Infinite.mono ?_ Nat.infinite_setOf_prime
  intro p hp
  obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A ⟨p, hp.pos⟩
  exact ⟨B, φ, prop_1_10_iv_mp P F hA φ ⟨hφ, by rw [hd]; exact hp⟩, by rw [hd]; rfl⟩

include P in
/-- ★★★**原文「there exist infinitely many isomorphism classes of objects of _A𝒞
that arise from irreducible arrows with domain A」そのもの**。

★★**検証役に 2 度差し戻された箇所である。**
1 度目は「『無限個』そのものを述べていない」、
2 度目は ★**「次数の集合の無限性(`_degrees_infinite`)と非同型性
(`_infinitely_many` の第2主張)は**別の集合**についての主張で、
橋渡しの一段がどの宣言にも無い」**。

★**私は docstring に「合わせて同型類の無限性にした」と書きながら、
その「合わせて」を書いていなかった。**
★★**「述べたと書いたものを述べていない」型の 4 例目**である。

★**ここで渡す**: `_A𝒞`(コスライス)の対象は `A` から出る射、同型は
「終域の同型 `θ` で `φ_m ≫ θ = φ_n` となるもの」。
素数の**単射列** `Nat.nth Nat.Prime` を次数に使えば、
次数が相異なるので `frobType_not_iso_of_degFr_ne` が非同型を与える。
★**`ℕ` で添字づけられた互いに非同型な族**が、「無限個の同型類」の忠実な形である。 -/
theorem prop_1_10_iv_iso_classes_infinite (F : FrobenioidCore P) {A : C}
    (hA : IsIsotropic P A) :
    ∃ f : ℕ → Σ B : C, A ⟶ B,
      (∀ n, IsIrreducibleMor (f n).2) ∧
      (∀ m n, m ≠ n → ∀ θ : (f m).1 ≅ (f n).1, (f m).2 ≫ θ.hom ≠ (f n).2) := by
  have hex : ∀ n : ℕ, ∃ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ ∧
      ((P.degFr φ : ℕ+) : ℕ) = Nat.nth Nat.Prime n := by
    intro n
    have hp : Nat.Prime (Nat.nth Nat.Prime n) := Nat.prime_nth_prime n
    obtain ⟨B, φ, hφ, hd⟩ := F.frobDegSurj A ⟨Nat.nth Nat.Prime n, hp.pos⟩
    exact ⟨B, φ, prop_1_10_iv_mp P F hA φ ⟨hφ, by rw [hd]; exact hp⟩, by rw [hd]; rfl⟩
  choose B φ hirr hdeg using hex
  refine ⟨fun n => ⟨B n, φ n⟩, hirr, fun m n hmn θ => ?_⟩
  haveI : IsIso θ.hom := inferInstance
  refine frobType_not_iso_of_degFr_ne P (φ m) (φ n) (fun h => hmn ?_) θ.hom
  refine Nat.nth_injective Nat.infinite_setOf_prime ?_
  rw [← hdeg m, ← hdeg n, h]

/-! ### ★(vi) 後半の最終段 —— `IsFrobeniusTrivial` を同型に沿って移す

★**原文(p.36)の最後の 1 文**:
> But by Proposition 1.4, (iii), these pre-steps are isomorphisms, so A is Frobenius-trivial.

★**「so A is Frobenius-trivial」の中身**は「同型に沿った移送」である。
★**mathlib に `End` の共役の準同型が無い**ので自作する
(`grep` で確認: `Mathlib/CategoryTheory/Endomorphism.lean` に `conj` は無い)。

★**`End` の乗法が `x * y = y ≫ x`** なので、共役 `f ↦ θ⁻¹ ≫ f ≫ θ` が
準同型であることの計算も**その向きで**行う。
-/

/-- ★**同型による `End` の共役**(`End A →* End B`)。

★mathlib に無いので置いた。★`End` の乗法は `x * y = y ≫ x` である。 -/
@[simps] def endConj {A B : C} (θ : A ≅ B) : End A →* End B where
  toFun f := θ.inv ≫ f ≫ θ.hom
  map_one' := by simp
  map_mul' x y := by
    simp only [End.mul_def]
    simp

include P in
/-- ★★**共役は `𝒪^×` を保つ**（2026-08-16 追加）。

★`Definition 1.2, (iv)` の `Frobenius-compact` の第3節で
`Aut_𝒞(C)` の `𝒪^×(C)` への作用を書くのに要る。

★**3 成分とも `Base` / `degFr` の関手性と同型の linear 性から出る**。 -/
theorem endConj_mem_otimes {A B : C} (θ : A ≅ B) {f : End A} (hf : f ∈ OTimes P A) :
    endConj θ f ∈ OTimes P B := by
  refine ⟨⟨?_, ?_⟩, IsUnit.map (endConj θ) hf.2⟩
  · show P.Base (θ.inv ≫ (f : A ⟶ A) ≫ θ.hom) = P.Base (𝟙 B)
    rw [P.Base_comp, P.Base_comp,
      show P.Base (f : A ⟶ A) = P.Base (𝟙 A) from hf.1.1,
      P.Base_id, Category.id_comp, ← P.Base_comp, θ.inv_hom_id]
  · show P.degFr (θ.inv ≫ (f : A ⟶ A) ≫ θ.hom) = 1
    rw [P.degFr_comp, P.degFr_comp,
      show P.degFr (f : A ⟶ A) = 1 from hf.1.2,
      show P.degFr θ.hom = 1 from isLinear_of_isIso P θ.hom,
      show P.degFr θ.inv = 1 from isLinear_of_isIso P θ.inv]
    simp

/-! ### ★★`Frobenius-compact`（2026-08-16 追加）

原文 (FrdI p.23):
> C such that O×(C) is commutative, O×(C)pf

★**引用を切った記録（事故 #3 の 7 度目）**: 続く `̸= 0, and every element of AutC(C)` は
★**`̸=`（打ち消し付きの等号）が抽出で落ちる**ため引用できない（35/63 文字で停止）。
★`′`・`▷`・`≼` に続く 4 つ目の「目では見えるが抽出に無い文字」である。

★★**`Pf` を経由せずに書く。** 理由は 2 つある。

1. `𝒪^×(C)` は **乗法**の `Submonoid (End A)` であり、`Pf` は
   `[AddCommMonoid M]` を要求する。★**しかも可換性は①として
   仮定される側**であって既知ではないので、型を作るのに仮定が要る。
2. ★★**同値な式がある**。

★**②の翻訳**（`Pf.mk_eq_zero_iff` が根拠）:
`M^pf = 0` ⇔ 「すべての `m` に `∃ k, k • m = 0`」。
乗法で書き直すと `M^pf ≠ 0` ⇔ ★**「ある `u ∈ 𝒪^×(C)` が無限位数を持つ」**。

★**③の翻訳**:
「`q = c/d` 倍として作用する」は `ℚ` を使わずに
`∀ x, d • σ x = c • x` と書ける。`M → M^pf` の像の上で書き直すと
`∃ k, (σ u)^(d k) = u^(c k)`。「acts trivially」は `∃ k, (σ u)^k = u^k`。

★★**像への制限が情報を落とさない理由**:
任意の `x = mk m a` は `a • x = of m` を満たす。誘導写像は加法的なので
`a • (d • σ x) = d • σ (a • x) = d • σ (of m) = c • (of m) = a • (c • x)`。
★**`Pf` では `a •` が単射**（`Pf.nsmul_injective`）なので、
像の上の式から `M^pf` 全体の式が出る。

★一度ここを「`Pf.eq_of_qact` が根拠」と書いたが、それは「同じ `c/d` 倍として
作用する写像は一意」という**別の主張**であり、ここでは効かない。
★★**検証役の再監査で訂正した** —— 根拠として挙げた補題と、実際に効く補題が違っていた。

★★**原文の三節をそのまま三連言にする。** -/

include P in
/-- **[FrdI] Definition 1.2, (iv)** `Frobenius-compact object`。

★上の docstring の翻訳に従う。★**共役が `𝒪^×` を保つ**こと
（`endConj_mem_otimes`）が③を書くのに要る。 -/
def IsFrobeniusCompact (A : C) : Prop :=
  (∀ x y : End A, x ∈ OTimes P A → y ∈ OTimes P A → x * y = y * x) ∧
  (∃ u : End A, u ∈ OTimes P A ∧ ∀ k : ℕ+, (u ^ ((k : ℕ+) : ℕ) : End A) ≠ 1) ∧
  (∀ (θ : A ≅ A) (c d : ℕ+),
    (∀ u : End A, u ∈ OTimes P A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) →
    ∀ u : End A, u ∈ OTimes P A → ∃ k : ℕ+,
      ((endConj θ u) ^ ((k : ℕ+) : ℕ) : End A) = (u ^ ((k : ℕ+) : ℕ) : End A))

include P in
/-- **[FrdI] Definition 1.2, (v)** `of Frobenius-compact type`。 -/
def IsOfFrobeniusCompactType : Prop := ∀ A : C, IsFrobeniusCompact P A

def IsFrobeniusCompact.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact object",
    sectionId := "frdi-def-1-2-iv" }

include P in
/-- ★★**`IsFrobeniusTrivial` は同型に沿って移る**。

★原文の「so A is Frobenius-trivial」の中身。 -/
theorem isFrobeniusTrivial_of_iso (F : FrobenioidCore P) {A B : C} (θ : A ≅ B)
    (h : IsFrobeniusTrivial P A) : IsFrobeniusTrivial P B := by
  obtain ⟨ζ, hd, hp⟩ := h
  refine ⟨(endConj θ).comp ζ, ?_, ?_⟩
  · intro n
    show P.degFr (θ.inv ≫ (ζ n : A ⟶ A) ≫ θ.hom) = n
    rw [P.degFr_comp, P.degFr_comp, degFr_of_isIso P θ.hom, degFr_of_isIso P θ.inv,
      hd n, one_mul, mul_one]
  · intro n
    constructor
    · show P.Base (θ.inv ≫ (ζ n : A ⟶ A) ≫ θ.hom) = P.Base (𝟙 B)
      rw [P.Base_comp, P.Base_comp, show P.Base (ζ n : A ⟶ A) = P.Base (𝟙 A) from (hp n).1,
        P.Base_id, Category.id_comp, ← P.Base_comp, θ.inv_hom_id, P.Base_id]
    · exact IsFrobeniusType.comp P F (isFrobeniusType_of_isIso P θ.inv)
        (IsFrobeniusType.comp P F (hp n).2 (isFrobeniusType_of_isIso P θ.hom))

include P in
/-- ★**group-like は base-isomorphism に沿って移る**。

★`Φ.map (Base α)` は `Base α` が同型なので全単射(`Φ.map (inv (Base α))` が逆)。
`Φ(A)` の元がすべて 0 なら `Φ(A′)` の元もすべて 0。 -/
theorem isGroupLikeObj_of_baseIso {A A' : C} (α : A' ⟶ A) (hbi : IsBaseIsomorphism P α)
    (hA : IsGroupLikeObj P A) : IsGroupLikeObj P A' := by
  haveI : IsIso (P.Base α) := hbi
  show IsGroupLike (Φ.val (P.toElem.obj A').base)
  rw [isGroupLike_iff]
  intro y
  have hy : y = Φ.map (P.Base α) (Φ.map (inv (P.Base α)) y) := by
    rw [← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
  have h0 : Φ.map (inv (P.Base α)) y = 0 :=
    eq_zero_of_isGroupLike_of_isSharp hA (P.divisorial _).2 _
  rw [hy, h0, map_zero]
  exact isAddUnit_zero

include P in
/-- ★★★**`Proposition 1.10, (vi)` 後半 完成** ——
isotropic 型で group-like な対象は Frobenius-trivial。

★原文(p.36)の最後の 3 文:
「`Definition 1.3, (i), (a), (b)` から co-angular pre-step `A′ → A`、`A′ → A″`
(`A″` は Frobenius-trivial)がある。`Proposition 1.4, (iii)` によりこれらは同型。
よって `A` は Frobenius-trivial。」

★**`A′` が group-like であること**(`isGroupLikeObj_of_baseIso`)が
★**原文が書いていない一歩**である —— pre-step が同型であることに要る isometry は
`A′` の group-like 性から来るのであって、`A` のそれからではない。 -/
theorem prop_1_10_vi_groupLike (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A A' A'' : C} (hA : IsGroupLikeObj P A)
    (α : A' ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (γ : A' ⟶ A'') (hγc : IsCoAngular P γ) (hγs : IsPreStep P γ)
    (hA'' : IsFrobeniusTrivial P A'') : IsFrobeniusTrivial P A := by
  -- `A′` も group-like
  have hA' : IsGroupLikeObj P A' := isGroupLikeObj_of_baseIso P α hαs.2 hA
  -- `α`・`γ` は同型
  haveI hαi : IsIso α := isIso_of_preStep_of_isGroupLikeObj' P F hiso hA' α hαs
  haveI hγi : IsIso γ := isIso_of_preStep_of_isGroupLikeObj' P F hiso hA' γ hγs
  -- `A″ ≅ A′ ≅ A`
  exact isFrobeniusTrivial_of_iso P F ((asIso γ).symm ≪≫ asIso α) hA''

/-! ### ★★(vi) の主語を置く —— `𝒞^istr` への具体化(2026-08-16)

原文 (FrdI p.35):
> (vi) The Frobenioid Cistr is of sub-quasi-Frobenius-trivial type. Moreover,

原文 (FrdI p.35):
> every group-like object A ∈Ob(Cistr) is Frobenius-trivial.

★★**取り下げ表の (vi) の残り 2 件**。`prop_1_10_vi_ofType` /
`prop_1_10_vi_groupLike` はどちらも「isotropic 型の Frobenioid」について述べており、
★**原文の主語 `𝒞^istr` がコメントの外に一度も現れていなかった。**

★**具体化そのものは 1 手である** —— `Proposition 1.9` が
`istr_frobenioid`(`𝒞^istr` は Frobenioid)と `istr_isotropic`(全対象が isotropic)を
既に与えているので、それを当てるだけでよい。
★★**「主語が不在」という欠落は、材料が揃っていても消えない** ——
一般形を証明したことと、原文が名指した対象について述べたことは別である。
-/

/-- ★★★**`Proposition 1.10, (vi)` 前半 —— 原文の主語で**。
`𝒞^istr` は sub-quasi-Frobenius-trivial 型である。 -/
theorem prop_1_10_vi_istr (F : FrobenioidCore P) (G : Frobenioid P) :
    IsOfSubQuasiFrobeniusTrivialType (istrPre P F) :=
  prop_1_10_vi_ofType (istrPre P F) (istr_frobenioid P F G) (istr_isotropic P F)

/-- ★★★**`Proposition 1.10, (vi)` 後半 —— 原文の主語で**。
`𝒞^istr` の group-like な対象は Frobenius-trivial である。

★`α`・`γ`・`Cc` は `prop_1_10_vi_data` が `Definition 1.3, (i), (a), (b)` から
導くので、ここでも仮定に置かない。 -/
theorem prop_1_10_vi_istr_groupLike (F : FrobenioidCore P) (G : Frobenioid P)
    (A : Istr P) (hA : IsGroupLikeObj (istrPre P F) A) :
    IsFrobeniusTrivial (istrPre P F) A := by
  obtain ⟨_, _, α, γ, hαc, hαs, hγc, hγs, hCc⟩ :=
    prop_1_10_vi_data (istrPre P F) (istr_frobenioid P F G) (istr_isotropic P F) A
  exact prop_1_10_vi_groupLike (istrPre P F) (istr_frobenioidCore P F)
    (istr_isotropic P F) hA α hαc hαs γ hγc hγs hCc

/-! ## ★★★出典の紐付け(`.src`) —— `Proposition 1.10` は **21 主張すべて完成**

★**`.src` は「その原典項目を完全に実装した」という主張である**(2026-08-15 に明文化した規則)。
`Proposition 1.10` は 6 条 21 主張:

| 条 | 主張 | 実装 |
|---|---|---|
| (i) | 10 | 一意性 / 存在(3 場合＋任意) / `degFr` / `Div` / 7 タイプ |
| (ii) | 3 | 組み替え / `degFr` / `Div` |
| (iii) | 3 | 像のモノイド / `𝒪^▷(A)` / `𝒪^×(A)` |
| (iv) | 2 | iff / 無限個の同型類 |
| (v) | 1 | iff |
| (vi) | 2 | sub-quasi-Frobenius-trivial / group-like ⟹ Frobenius-trivial |

★★**この数え方が誤っていた**(下の取り下げ記録を見よ)。
-/

/-! ## ★★`.src` を**取り下げた**記録(2026-08-16)

★私の文脈を持たない検証役(Challenger)に監査を依頼し、
★★**6 個すべてが規則(「`.src` は完全実装の主張」)を満たしていない**ことが判明した。
★指摘はすべて**原文 (p.34–35) と Lean を私自身が直接読んで確認済み**である。

| 条 | 欠けているもの | 該当宣言 |
|---|---|---|
| (i) | ~~★`β` の量化子が逆~~ → ★**実装した**。`Definition 1.3, (ii)` の本質的一意性（`frobDegUniq`）で、自分で作った `β₀` を与えられた `β` に合わせる | `prop_1_10_i_exists_given` |
| (i) | ~~原文「In this situation, degFr(φ) = degFr(φ′)」が**ファイルに存在しない**~~ → ★**実装した** | `prop_1_10_i_degFr_phi_eq` |
| (i) | ★★原文「then the same is true of **φ′**」の 7 タイプ。`prop_1_10_i_four_types` は `φ` についての主張で `φ′` のものではない。★★**4 つすべて実装した**(`prop_1_10_i_baseIso_of` / `_isometric_of` / `_coAngular_of` / `_lbInvertible_of`)。★co-angular の鍵は「`φ ≫ β` が co-angular」であることだった —— 引き戻す必要はなく、**分解の側を前合成で延ばせばよかった**。★★**7 タイプすべて実装完了**（上の 4 本 ＋ `prop_1_10_i_linear_of` / `prop_1_10_i_preStep_of` / `prop_1_10_i_frobType_of` / `prop_1_10_i_pullBack_of`）。★pull-back は `Proposition 1.4, (ii)` を両向きに使って普遍性を避けた | `prop_1_10_i_four_types` |
| (ii) | ~~`Div` の式が `β′` の base-isomorphism 性を仮定せず、原文の `β′∗`(全単射)の形になっていない~~ → ★**実装した**。`β′` は Frobenius 型ゆえ base-isomorphism なので、`Φ.map (Base β′)` の逆で原文の向きに直せる | `prop_1_10_ii_Div_formula'` |
| (iii) | ~~★原文は証明中で「`A` を Frobenius-trivial としてよい」と**還元している**が、その還元が未実装で、仮定に逃がしている~~ → ★**実装した**。原文が名指す 3 つ(`Definition 1.3, (i), (a)` / `(i), (b)` / `(iii), (c)`)がそのまま 3 段になった。★`isotropic` の最初の用途が**span の pre-step を co-angular にすること**だったと分かった | `prop_1_10_iii_otri_perfect_of_type`, `prop_1_10_iii_otimes_perfect_of_type` |
| (iii) | ~~「the monoids in the image of Φ」は Ob(𝒟) 全体の像だが、実装は `Base` の像のみ~~ → ★**実装した**。`baseSurj` の同型を `MonoidOn.map_bijective_of_iso` で渡り、perfect 性を `isPerfectMonoid_of_bijective` で移した | `prop_1_10_iii_image_perfect_of_base` |
| (iv) | ~~★原文は「**域が** isotropic」だけを要求するが、実装は `∀ X : C, IsIsotropic P X`(**圏全体**)を仮定~~ → ★**実装した**。`isotropicClosed` で域の isotropic 性から下流へ伝播させれば足りた | `prop_1_10_iv`, `prop_1_10_iv_mp` |
| (iv) | ~~「infinitely many isomorphism classes」そのものを述べていない(素数ごとの存在と非同型性まで)~~ → ★**実装した**。★**`ℕ` で添字づけられた互いに非同型な既約射の族**として述べた。★★**当初は「次数の集合が無限」＋「次数が違えば非同型」で済ませ、表にも「合わせて同型類の無限性にした」と書いたが、その「合わせて」の一手をどの宣言にも書いていなかった**(検証役が発見) | `prop_1_10_iv_iso_classes_infinite` |
| (v) | ~~`IsPrimeFrobComposite` の基底を「任意の同型」に取った——**意図的な修復**だが、原文項目そのものの実装ではない~~ → ★★**反例に格上げした**。素直な読み(基底＝恒等射)を `IsPrimeFrobCompositeId` として別に定義し、★**恒等射でない同型はその形に書けない**ことを示した。同型は Frobenius 型なので、素直な読みでは iff が破れる —— ★**基底を同型に取ったのは修復ではなく強制であった** | `not_isPrimeFrobCompositeId_of_isIso_of_ne_id` |
| (vi) | ~~★★**`𝒞^istr` がコメント外に一度も現れない**。(vi) の主語が不在~~ → ★**実装した**。`Proposition 1.9` の `istr_frobenioid` / `istr_isotropic` を当てるだけだった | `prop_1_10_vi_istr` |
| (vi) | ~~★`α`,`γ`,`Cc` の存在を**仮定に置いている**が、原文は `Definition 1.3, (i), (a), (b)` から**導いている**~~ → ★**実装した** | `prop_1_10_vi_data`, `prop_1_10_vi_ofType` |
| (vi) | ~~「Moreover, every group-like object … is Frobenius-trivial.」が別宣言で、そちらも `𝒞^istr` への具体化がない~~ → ★**実装した** | `prop_1_10_vi_istr_groupLike` |

★★**ここにある定理はどれも真であり、`sorry` も無い。**
★失われたのは「完全である」という**主張**だけである。
★**上の表がそのまま次の作業のチェックリストになった。**

★★**2026-08-16: 表の 12 行すべてが埋まった。**
★ただし ★**`.src` を付けるのは、私の文脈を持たない検証役が改めて全条を監査してからである** ——
過去 3 回、「埋めた」と数えた項目に**別種の**欠落が見つかっている。
★**表が埋まったことは「もう欠落が無い」ことを意味しない。**

★★**教訓**: 自分の文脈を継いだ検証は、自分の誤読も継ぐ。
「21/21」と数えたのは私であり、その数え方は「**実装したものを数えた**」ものであって、
★**「原文が述べたものを数えた」ものではなかった。**

## ★★★2026-08-16 —— 取り下げから復帰した

★**私の文脈を持たない検証役が、原文 (p.34–35) を自分で読んで主張を数え直し
(22 主張。`∃!` を存在と一意性に分けたため。1 と数えれば 21 で一致する)、
22 すべてに対応する宣言があることを確認した。**

★**過去に落ちた 5 つの型を個別に検査して、残っているものは無い**:
量化子の向き / 原文より強い仮定 / 省いた段を仮定に逃がす / 主語の不在 /
docstring の理由違い。

★★**監査は 2 往復した。** 1 往復目の判定は「いけない」で、理由は
★**「次数の集合の無限性」と「次数が違えば非同型」は別の集合についての主張で、
橋渡しの一段がどの宣言にも無い」**だった ——
★**私は表に「合わせて同型類の無限性にした」と書きながら、その「合わせて」を
書いていなかった。**「述べたと書いたものを述べていない」型の 4 例目である。

★**ここで条なしの `.src` を付ける。**
-/

/-- ★★★**[FrdI] Proposition 1.10 全体**。(i)〜(vi) の 22 主張がすべて実装された。 -/
def prop_1_10.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 34, item := "Proposition 1.10",
    sectionId := "frdi-prop-1-10" }

/-! ## ★★`base-trivial` な対象は `Frobenius-trivial`(2026-08-16 追加)

★★**`Proposition 1.6, (v)` の `⟸` を証明する道の第 1 歩**である。
`Gap_1_6_v` に不足しているのは `Aut-ample`(底を指定した同型の取り直し)だが、
★**反例の設計が 2 度とも `Definition 1.3` に弾かれた**ので、
★**「出ない」ではなく「出る」方を試す。**

★中身は `Definition 1.3, (i), (a)` そのもの ——
`𝒟` のどの同型類も **Frobenius-trivial な対象**から来るので、
`base-trivial` はその対象を自分自身へ引き寄せる。

★**道具は既にあった**(`isFrobeniusTrivial_of_iso`)—— 2026-08-16 に
「無い」と思って書き直しかけたが、`grep` の出力を `head` で切っていただけだった。 -/

include P in
/-- ★★★**`base-trivial` な対象は `Frobenius-trivial`**。

★`Definition 1.3, (i), (a)` が「底の同型類は Frobenius-trivial な対象から来る」と言い、
`base-trivial` が「底が同型なら対象も同型」と言う。★**噛み合う。** -/
theorem isFrobeniusTrivial_of_baseTrivial (F : FrobenioidCore P) {A : C}
    (h : IsBaseTrivial P A) : IsFrobeniusTrivial P A := by
  obtain ⟨A₀, hft, ⟨e⟩⟩ := F.baseSurj (P.toElem.obj A).base
  obtain ⟨k⟩ := h A₀ ⟨e.symm⟩
  exact isFrobeniusTrivial_of_iso P F k hft

include P in
/-- ★★★**pre-step 自己射は、`Base` を保ったまま `Div` を `n` 倍できる**。

★★**`Proposition 1.6, (v)` の `⟸` を証明する道の第 3 歩**である。

★**中身**: `A` が Frobenius-trivial なら底恒等な Frobenius 自己射 `ζ_n` がある。
`s ≫ ζ_n` を `Definition 1.3, (iv), (a)` で 3 分解し、Frobenius 因子を
(ii) の本質的一意性で `ζ_n` に取り替えると `s ≫ ζ_n = ζ_n ≫ t` になる。
★**`ζ_n` は底恒等で `Div = 0` なので、`Base t = Base s` かつ `Div t = n · Div s`。** -/
theorem preStep_endo_scale (F : FrobenioidCore P) {A : C} (hft : IsFrobeniusTrivial P A)
    (s : A ⟶ A) (hs : IsPreStep P s) (n : ℕ+) :
    ∃ t : A ⟶ A, IsPreStep P t ∧ P.Base t = P.Base s ∧
      P.Div t = ((n : ℕ+) : ℕ) • P.Div s := by
  obtain ⟨ζ, hdeg, hprop⟩ := hft
  have hzb : P.Base (ζ n : A ⟶ A) = 𝟙 _ := by
    have h := (hprop n).1
    show P.Base (ζ n : A ⟶ A) = 𝟙 _
    rw [← P.Base_id A]
    exact h
  have hzd : P.Div (ζ n : A ⟶ A) = 0 := (hprop n).2.1.2
  -- ★`s ≫ ζ n` を `Definition 1.3, (iv), (a)` で 3 分解する
  obtain ⟨X, Y, γ, β, α, hx, hγ, hβ, hα⟩ := F.arbFactor (s ≫ (ζ n : A ⟶ A))
  obtain ⟨-, hαl⟩ := F.pullBackLB α hα
  have hdγ : P.degFr γ = n := by
    have h := congrArg P.degFr hx
    rw [P.degFr_comp, hdeg n, hs.1, mul_one, P.degFr_comp, P.degFr_comp, hαl, hβ.1,
      one_mul, one_mul] at h
    exact h.symm
  -- ★(ii) の本質的一意性で Frobenius 因子を `ζ n` に取り替える
  obtain ⟨θ, hθi, hθ⟩ := F.frobDegUniq A X A γ (ζ n : A ⟶ A) hγ (hprop n).2
    (by rw [hdγ, hdeg n])
  haveI := hθi
  refine ⟨inv θ ≫ β ≫ α, ⟨?_, ?_⟩, ?_, ?_⟩ <;>
    [skip; skip; skip; skip]
  case refine_1 =>
    show P.degFr (inv θ ≫ β ≫ α) = 1
    rw [P.degFr_comp, P.degFr_comp, hαl, hβ.1, degFr_of_isIso P (inv θ)]
    simp
  all_goals
    have hxt : s ≫ (ζ n : A ⟶ A) = (ζ n : A ⟶ A) ≫ (inv θ ≫ β ≫ α) := by
      rw [hx, ← hθ]; simp
    have hbase : P.Base (inv θ ≫ β ≫ α) = P.Base s := by
      have h := congrArg P.Base hxt
      rw [P.Base_comp, P.Base_comp, hzb, Category.comp_id, Category.id_comp] at h
      exact h.symm
  case refine_2 =>
    show IsIso (P.Base (inv θ ≫ β ≫ α))
    rw [hbase]
    exact hs.2
  case refine_3 => exact hbase
  case refine_4 =>
    have h := congrArg P.Div hxt
    rw [P.Div_comp, P.Div_comp, hzb, hzd, hdeg n] at h
    simp only [map_zero, zero_add, smul_zero, add_zero, MonoidOn.map_id, id_eq] at h
    exact h.symm

include P in
/-- ★★**自己射はすべて co-angular** —— `Definition 1.3, (iii), (b)` を `𝟙` に当てるだけ。

★`𝟙_A` は co-angular pre-step なので、(iii)(b) が `A ⟶ A` のすべてを co-angular にする。 -/
theorem endo_isCoAngular (F : FrobenioidCore P) {A : C} (φ : A ⟶ A) : IsCoAngular P φ :=
  F.coAngularOfPreStep (𝟙 A) (isCoAngular_id P A) (isPreStep_id P A) φ

include P in
/-- ★★★**`base-trivial` ⟹ `Aut-ample`**(★捻れが自明な場合)。

★★**`Proposition 1.6, (v)` の `⟸` を証明する道の組み立て**である。

★**筋**:
1. `preStepSpan` が `α = (Base φ)⁻¹ ≫ Base ψ`(`φ, ψ` は pre-step)を与える
2. `base-trivial` で中間対象を `A` へ寄せ、`φ, ψ` を **`A` の自己射**にする
3. 自己射はすべて co-angular(`endo_isCoAngular`)
4. `preStep_endo_scale` で `Div` を 2 倍した `t`(`Base t = Base s`)を作り、
   `s ≫ s`(`Div = Φ(Base s)(Div s) + Div s`)と `Div` を突き合わせる
5. `coaPre_base_diff` が同型 `f` を与え、`Base f = Base s` が出る

★★**残る仮定 `htwist` が唯一の穴である** ——
`Φ(Base s)(Div s) = Div s`、すなわち ★**pre-step 自己射の底が自分の `Div` を動かさない**こと。
`Φ` が定数関手ならこれは自明だが、一般には言えていない。
★**これが `Gap_1_6_v` の正体を 1 点に絞ったものである。** -/
theorem isAutAmple_of_baseTrivial_of_untwisted [MorphismProperty.IsMultiplicative (coaPreProp P)]
    (F : FrobenioidCore P)
    (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence)
    {A : C} (h : IsBaseTrivial P A)
    (htwist : ∀ s : A ⟶ A, IsPreStep P s → Φ.map (P.Base s) (P.Div s) = P.Div s) :
    IsAutAmple P A := by
  have hft : IsFrobeniusTrivial P A := isFrobeniusTrivial_of_baseTrivial P F h
  -- ★核心: pre-step 自己射の底は、同型の底として実現できる
  have key : ∀ s : A ⟶ A, IsPreStep P s → ∃ u : A ⟶ A, IsIso u ∧ P.Base u = P.Base s := by
    intro s hs
    obtain ⟨t, hts, htb, htd⟩ := preStep_endo_scale P F hft s hs 2
    have hdiv : P.Div t = P.Div (s ≫ s) := by
      rw [htd, P.Div_comp, htwist s hs, hs.1]
      simp [two_nsmul]
    obtain ⟨f, hfi, hf⟩ := coaPre_base_diff P hequiv t (s ≫ s)
      ⟨endo_isCoAngular P F t, hts⟩ ⟨endo_isCoAngular P F _, hs.comp P hs⟩ hdiv
    haveI := hfi
    haveI : IsIso (P.Base s) := hs.2
    refine ⟨f, hfi, ?_⟩
    have hb := congrArg P.Base hf
    rw [P.Base_comp, P.Base_comp, htb] at hb
    exact (cancel_epi (P.Base s)).mp hb
  intro g hg
  haveI := hg
  obtain ⟨X, φ, ψ, hφ, hψ, hgeq⟩ := F.preStepSpan A A g hg
  haveI hφb : IsIso (P.Base φ) := hφ.2
  obtain ⟨e⟩ := h X ⟨(asIso (P.Base φ)).symm⟩
  -- ★`e : X ≅ A`。`e.inv ≫ φ`、`e.inv ≫ ψ` は `A` の pre-step 自己射
  have hφ' : IsPreStep P (e.inv ≫ φ) := (isPreStep_of_isIso P e.inv).comp P hφ
  have hψ' : IsPreStep P (e.inv ≫ ψ) := (isPreStep_of_isIso P e.inv).comp P hψ
  obtain ⟨u, hui, hub⟩ := key _ hφ'
  obtain ⟨v, hvi, hvb⟩ := key _ hψ'
  haveI := hui
  haveI := hvi
  refine ⟨inv u ≫ v, inferInstance, ?_⟩
  haveI hei : IsIso (P.Base e.inv) := isBaseIsomorphism_of_isIso P e.inv
  haveI hbui : IsIso (P.Base u) := isBaseIsomorphism_of_isIso P u
  have hbu : P.Base u = P.Base e.inv ≫ P.Base φ := by rw [hub, P.Base_comp]
  have hbv : P.Base v = P.Base e.inv ≫ P.Base ψ := by rw [hvb, P.Base_comp]
  -- ★`inv` の下を書き換えると motive が壊れるので、`P.Base u` を前から掛けて消す
  have hpu : P.Base u ≫ P.Base (inv u ≫ v) = P.Base v := by
    rw [← P.Base_comp, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  refine (cancel_epi (P.Base u)).mp ?_
  rw [hpu, hbv, hbu, hgeq]
  simp

include P in
/-- ★★★**`base-trivial` ⟹ `Aut-ample`**(★捻れに触れない版)。

★★**`htwist` を使わない** —— 代わりに
`hsurj`「`Φ(A)` のどの元も **底恒等な pre-step 自己射の `Div`** として現れる」
だけを仮定する。

★★**筋がずっと素直になる**:
`preStepSpan` の 2 本 `φ', ψ'`(`A` の pre-step 自己射)の `Div` を
**後から底恒等な pre-step を継いで揃える** ——
`Φ(Base φ')` は全単射(`Definition 1.1, (ii), (b)`)なので
`Div r₁ := Φ(Base φ')⁻¹(Div ψ')` が取れて、
`Div (φ' ≫ r₁) = Div ψ' + Div φ'`。同様に `Div (ψ' ≫ r₂) = Div φ' + Div ψ'`。
★**底は `𝟙` を継ぐので変わらず、比 `g` も保たれる。**
あとは `coaPre_base_diff` が同型 `f` を与え、★**`Base f = g` がそのまま出る。**

★★**正直な注記**: `hsurj` は `base-trivial` の下で **`Aut-ample` と同値**である
(逆向きは、`Div = x` の自己射の底を `Aut-ample` で打ち消せばよい)。
★**したがってこれは穴を閉じたのではなく、穴を言い換えたものである。**
★ただし `htwist`(捻れの技術的条件)より **`Definition 1.2` の語彙に近い**形になった:
「**`𝒪^▷(A) → Φ(A)` が全射**」。 -/
theorem isAutAmple_of_baseTrivial_of_divSurj
    [MorphismProperty.IsMultiplicative (coaPreProp P)] (F : FrobenioidCore P)
    (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence)
    {A : C} (h : IsBaseTrivial P A)
    (hsurj : ∀ x : Φ.val (P.toElem.obj A).base, ∃ r : A ⟶ A, IsPreStep P r ∧
      P.Base r = 𝟙 _ ∧ P.Div r = x) :
    IsAutAmple P A := by
  intro g hg
  haveI := hg
  obtain ⟨X, φ, ψ, hφ, hψ, hgeq⟩ := F.preStepSpan A A g hg
  haveI hφb : IsIso (P.Base φ) := hφ.2
  obtain ⟨e⟩ := h X ⟨(asIso (P.Base φ)).symm⟩
  set φ' : A ⟶ A := e.inv ≫ φ with hφ'def
  set ψ' : A ⟶ A := e.inv ≫ ψ with hψ'def
  have hφ' : IsPreStep P φ' := (isPreStep_of_isIso P e.inv).comp P hφ
  have hψ' : IsPreStep P ψ' := (isPreStep_of_isIso P e.inv).comp P hψ
  haveI hb1 : IsIso (P.Base φ') := hφ'.2
  haveI hb2 : IsIso (P.Base ψ') := hψ'.2
  -- ★`Φ(Base φ')` は全単射
  obtain ⟨y₁, hy₁⟩ := (Φ.fsmIso (P.Base φ') (isFSMMorphism_of_isIso _)).2 (P.Div ψ')
  obtain ⟨y₂, hy₂⟩ := (Φ.fsmIso (P.Base ψ') (isFSMMorphism_of_isIso _)).2 (P.Div φ')
  have hy₁' : Φ.map (P.Base φ') y₁ = P.Div ψ' := hy₁
  have hy₂' : Φ.map (P.Base ψ') y₂ = P.Div φ' := hy₂
  obtain ⟨r₁, hr₁s, hr₁b, hr₁d⟩ := hsurj y₁
  obtain ⟨r₂, hr₂s, hr₂b, hr₂d⟩ := hsurj y₂
  -- ★`Div` を揃える
  have hd1 : P.Div (φ' ≫ r₁) = P.Div ψ' + P.Div φ' := by
    rw [P.Div_comp, hr₁d, hy₁', hr₁s.1]
    simp
  have hd2 : P.Div (ψ' ≫ r₂) = P.Div φ' + P.Div ψ' := by
    rw [P.Div_comp, hr₂d, hy₂', hr₂s.1]
    simp
  have hdeq : P.Div (φ' ≫ r₁) = P.Div (ψ' ≫ r₂) := by
    rw [hd1, hd2, add_comm]
  -- ★底は変わらない
  have hb1' : P.Base (φ' ≫ r₁) = P.Base φ' := by rw [P.Base_comp, hr₁b, Category.comp_id]
  have hb2' : P.Base (ψ' ≫ r₂) = P.Base ψ' := by rw [P.Base_comp, hr₂b, Category.comp_id]
  obtain ⟨f, hfi, hf⟩ := coaPre_base_diff P hequiv (φ' ≫ r₁) (ψ' ≫ r₂)
    ⟨endo_isCoAngular P F _, hφ'.comp P hr₁s⟩ ⟨endo_isCoAngular P F _, hψ'.comp P hr₂s⟩ hdeq
  refine ⟨f, hfi, ?_⟩
  have hbf := congrArg P.Base hf
  rw [P.Base_comp, hb1', hb2'] at hbf
  -- ★`Base φ' ≫ Base f = Base ψ'`、かつ `g = (Base φ')⁻¹ ≫ Base ψ'`
  have hbφ' : P.Base φ' = P.Base e.inv ≫ P.Base φ := by rw [hφ'def, P.Base_comp]
  have hbψ' : P.Base ψ' = P.Base e.inv ≫ P.Base ψ := by rw [hψ'def, P.Base_comp]
  haveI hei : IsIso (P.Base e.inv) := isBaseIsomorphism_of_isIso P e.inv
  refine (cancel_epi (P.Base φ')).mp ?_
  rw [hbf, hgeq, hbφ', hbψ']
  simp

end ABC3.Found.FrdI
