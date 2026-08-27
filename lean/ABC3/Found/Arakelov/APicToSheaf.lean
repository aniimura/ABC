/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.SheafifyTriv
import ABC3.Found.Arakelov.PicEvalIso
import ABC3.Found.Arakelov.PicSheafGroup
import ABC3.Found.Arakelov.PicClassGroup
import ABC3.Found.Arakelov.PicAssoc
import ABC3.Found.Arakelov.PicGammaInv
import Mathlib.Algebra.Module.PID
import ABC3.Meta.Claim

/-!
# **`APic`（前層・計量つき）から `Pic`（層）へ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★台帳 `arakelov-degF-finite-places` の**段 A3**

無条件の `deg_F`（[Szp] を仮定しない形）に要るのは

    `ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)`

だけになった（`§9-746`・`§9-747`）。その**有限素点側**は

    `APicM (Spec 𝓞_F) → PicSheaf (Spec 𝓞_F) ≃* Cl(F)`

で読める。★★本ファイルはその**第 1 の矢印**を作る。

| 段 | 内容 | 状態 |
|---|---|---|
| A1 | 層化は局所自明性を保つ | ✅ 在庫 |
| A2 | 自明化・遷移単元の輸送 | ✅ `SheafifyTriv.lean`（§9-748） |
| A3 | **`APicM X →* PicSheaf X`** | ★★**本ファイル**（群準同型） |
| — | `PicSheaf (Spec 𝓞_F) ≃* Cl(F)` | ✅ 在庫（`picSheafEquivClassGroupOF`） |

## ★★★★★★★★在庫が全部持っていた

★測定（2026-08-28）で分かったのは、要る部品が**すべて在庫**だったことである:

| 部品 | 場所 |
|---|---|
| `isLocallyTrivial_sheafify`（層化は局所自明性を保つ） | `PicSheafifyTrivial.lean` |
| `InvSheaf.ofLocallyTrivial`（局所自明な層加群は可逆層） | `PicEvalIso.lean` |
| `PicSheaf.mk_eq_mk`（同型なら同じ類） | `PicSheafGroup.lean` |
| `picSheafEquivClassGroupOF`（`Pic(Spec 𝓞_F) ≃* Cl(F)`） | `PicClassGroup.lean` |

★★したがって本ファイルは**それらを継ぐだけ**である。

## ★★★★★★★★群準同型であること —— 層化はテンソル積と可換

`map_mul` に要るのは

    `層化 (A ⊗ B) ≅ tensorModules (層化 A) (層化 B)`

であり、★これも在庫を 2 本継ぐだけであった:

    `sheafifyTensorRight`：`層化 (A ⊗ B) ≅ 層化 ((層化 A).val ⊗ B)`
    `sheafifyTensorLeft` ：`層化 (A' ⊗ B) ≅ 層化 (A' ⊗ (層化 B).val)`

★★`map_one` は要らない——`APicM X` も `PicSheaf X` も**群**なので
`MonoidHom.mk'` が `map_mul` だけから作る。

## ★残っている段（明示）

★★アルキメデス成分（計量から無限素点の実数を読む段）と、
`ADiv/APrc ≅ APicM` の組み立て（逆向きを含む）が残る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★層化して可逆層と見る -/

/-- ★★★★★★**計量つき算術直線束を層の水準の可逆層と見る**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★層化して `InvSheaf.ofLocallyTrivial` に渡すだけである
——局所自明性は `isLocallyTrivial_sheafify`（在庫）が保つ。 -/
noncomputable def AInv.toInvSheaf {X : Scheme.{0}} (L : AInv X) : InvSheaf X :=
  InvSheaf.ofLocallyTrivial ((sheafifyFunctor X).obj L.carrier.sheaf)
    (isLocallyTrivial_sheafify X L.carrier.sheaf L.carrier.triv)

/-- ★★★**等長同型なら層化して同じ類になる**。

★等長性は捨てて同型だけを使う——`Pic` は計量を見ないからである。 -/
theorem picSheaf_mk_eq_of_isometric {X : Scheme.{0}} {L M : AInv X}
    (h : Isometric L.carrier M.carrier) :
    PicSheaf.mk (AInv.toInvSheaf L) = PicSheaf.mk (AInv.toInvSheaf M) :=
  h.elim fun φ _ => PicSheaf.mk_eq_mk ((sheafifyFunctor X).mapIso φ)

/-! ## ★★★★★★★★`APicM X → PicSheaf X` -/

/-- ★★★★★★★★**`APic`（前層・計量つき）から `Pic`（層）への写像**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これが台帳 `arakelov-degF-finite-places` の**段 A3**（写像の水準）である。
★★★群準同型版は下の `APicMToPicSheafHom`。 -/
noncomputable def APicMToPicSheaf {X : Scheme.{0}} : APicM X → PicSheaf X :=
  Quotient.lift (fun L => PicSheaf.mk (AInv.toInvSheaf L))
    (fun _ _ h => picSheaf_mk_eq_of_isometric h)

@[simp] theorem APicMToPicSheaf_mk {X : Scheme.{0}} (L : AInv X) :
    APicMToPicSheaf (APicM.mk L) = PicSheaf.mk (AInv.toInvSheaf L) := rfl

/-! ## ★★★★★★★★★★有限素点側——類群へ -/

/-- ★★★★★★★★★★**`APicM (Spec 𝓞_F) → Cl(F)`**——`deg_F` の橋の有限素点側。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★在庫の `picSheafEquivClassGroupOF`（`Pic(Spec 𝓞_F) ≃* Cl(F)`）を継ぐだけである。
★★群準同型版は下の `APicMToClassGroupHom`。 -/
noncomputable def APicMToClassGroup (F : Type) [Field F] [NumberField F] :
    APicM (Spec (CommRingCat.of (𝓞 F))) → ClassGroup (𝓞 F) :=
  fun L => picSheafEquivClassGroupOF F (APicMToPicSheaf L)

@[simp] theorem APicMToClassGroup_mk (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    APicMToClassGroup F (APicM.mk L)
      = picSheafEquivClassGroupOF F (PicSheaf.mk (AInv.toInvSheaf L)) := rfl

/-! ## ★★★★★★★★層化はテンソル積と可換である -/

/-- ★★★★★★★★**層化はテンソル積と可換である**（可逆層について）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★在庫の `sheafifyTensorRight` と `sheafifyTensorLeft` を継ぐだけである。 -/
noncomputable def sheafifyTensorIso {X : Scheme.{0}} (A B : X.PresheafOfModules)
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) :
    (sheafifyFunctor X).obj (A ⊗ B)
      ≅ tensorModules ((sheafifyFunctor X).obj A) ((sheafifyFunctor X).obj B) :=
  sheafifyTensorRight X A B hB.isLocallyRankOne
    ≪≫ sheafifyTensorLeft X ((sheafifyFunctor X).obj A).val B
        (isLocallyTrivial_sheafify X A hA).isLocallyRankOne

theorem APicMToPicSheaf_mul_mk {X : Scheme.{0}} (L M : AInv X) :
    APicMToPicSheaf (APicM.mk L * APicM.mk M)
      = APicMToPicSheaf (APicM.mk L) * APicMToPicSheaf (APicM.mk M) := by
  show PicSheaf.mk (AInv.toInvSheaf (L.mul M))
    = PicSheaf.mk ((AInv.toInvSheaf L).mul (AInv.toInvSheaf M))
  exact PicSheaf.mk_eq_mk
    (sheafifyTensorIso L.carrier.sheaf M.carrier.sheaf L.carrier.triv M.carrier.triv)

/-- ★★★★★★★★★★**`APicM X →* PicSheaf X` は群準同型である**——台帳の段 A3。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`map_one` は要らない——どちらも**群**なので `MonoidHom.mk'` が
`map_mul` だけから作る。 -/
noncomputable def APicMToPicSheafHom {X : Scheme.{0}} : APicM X →* PicSheaf X :=
  MonoidHom.mk' APicMToPicSheaf (by
    refine Quotient.ind fun L => Quotient.ind fun M => ?_
    exact APicMToPicSheaf_mul_mk L M)

/-- ★★★★★★★★★★**`APicM (Spec 𝓞_F) →* Cl(F)`**——`deg_F` の橋の有限素点側、
**群準同型として**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これで `ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)` の**有限素点側が群準同型として存在する**。
★★★残るのはアルキメデス成分と、逆向きの構成である。 -/
noncomputable def APicMToClassGroupHom (F : Type) [Field F] [NumberField F] :
    APicM (Spec (CommRingCat.of (𝓞 F))) →* ClassGroup (𝓞 F) :=
  (picSheafEquivClassGroupOF F).toMonoidHom.comp APicMToPicSheafHom

@[simp] theorem APicMToClassGroupHom_mk (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    APicMToClassGroupHom F (APicM.mk L)
      = picSheafEquivClassGroupOF F (PicSheaf.mk (AInv.toInvSheaf L)) := rfl

/-! ## ★★★★★★★★★★段 A の到達点 —— `Γ(L)` は可逆加群である -/

/-- ★★★★★★★★★★**台帳の段 A（アフィンの橋）そのもの**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

台帳 `arakelov-degF-finite-places` の段 A の claim は

> `Spec R` 上の局所自明な**前層**加群 `L` から**可逆 `R`-加群** `Γ(L,⊤)` を得る

であった。★層化して（`AInv.toInvSheaf`）在庫の `invertible_gammaCarrier` に渡すだけである。

★★これで**段 A が閉じた**。★★★段 B（有限部分 `log #(Γ(L)/𝓞_F·s)`）の前提が揃う
——`Γ(L)` が可逆（＝階数 1 射影）なので、`s ≠ 0` に対して商が有限になる。 -/
theorem invertible_gamma_toInvSheaf (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Invertible (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) :=
  invertible_gammaCarrier R (AInv.toInvSheaf L)

/-! ## ★★★★★段 B への一歩 —— 非零切断の存在 -/

open scoped TensorProduct in
/-- ★★★**可逆加群は自明でない**。

★`Mv ⊗ M ≃ R` なので、`M = 0` なら `R = 0` になる。 -/
theorem invertible_nontrivial (R M : Type) [CommRing R] [AddCommGroup M] [Module R M]
    [Module.Invertible R M] [Nontrivial R] : Nontrivial M := by
  by_contra h
  rw [not_nontrivial_iff_subsingleton] at h
  haveI : Subsingleton (Module.Dual R M ⊗[R] M) := by infer_instance
  haveI : Subsingleton R := (Module.Invertible.linearEquiv R M).symm.toEquiv.subsingleton
  exact false_of_nontrivial_of_subsingleton R

/-- ★★★★★**算術直線束には非零の大域切断がある**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★古典的な `deg_F` の式は **非零の `s`** を取る。
★★段 A が閉じたのでその存在が言える。 -/
theorem exists_ne_zero_gamma (R : CommRingCat.{0}) [Nontrivial (R : Type)]
    (L : AInv (Spec R)) :
    ∃ s : (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type),
      s ≠ 0 := by
  haveI := invertible_gamma_toInvSheaf R L
  haveI := invertible_nontrivial (R : Type)
    (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type)
  exact exists_ne 0

/-! ## ★★★★段 B の前提 —— 有限生成と射影性 -/

/-- ★★★**`Γ(L)` は有限生成である**。

★mathlib の「An invertible module is finite and projective」が
**インスタンスとして**与えてくれる(実測 2026-08-28)。 -/
theorem finite_gamma (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Finite (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) := by
  haveI := invertible_gamma_toInvSheaf R L
  infer_instance

/-- ★★★**`Γ(L)` は射影である**。★射影＋整域なので**捻れなし**であり、
これが段 B(商の有限性)の出発点になる。 -/
theorem projective_gamma (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Projective (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) := by
  haveI := invertible_gamma_toInvSheaf R L
  infer_instance

/-! ## ★★★★★★★段 B の核へ —— 分数体上で rank 1 -/

open scoped TensorProduct in
/-- ★★★★★★★**分数体に上げると rank 1 になる**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★実測(2026-08-28)で道が開けた——mathlib の
`Module.Invertible.finrank_eq_one` は `[Free R M]` を前提とするので
Dedekind 域上では直接使えないが、

  ・`instance : Module.Invertible A (A ⊗[R] M)`(環を変える底変換)と
  ・体上の加群は常に free

がどちらも**インスタンスで動く**ので、分数体へ上げれば 1 行で出る。

★★★★これが段 B(商の有限性)の鍵である
——rank 1 なので `s ≠ 0` は分数体上で全体を張り、
商 `Γ(L)/(R·s)` が捻れになる。 -/
theorem finrank_fractionRing_gamma (R : CommRingCat.{0}) [IsDomain (R : Type)]
    (L : AInv (Spec R)) :
    Module.finrank (FractionRing (R : Type))
        (FractionRing (R : Type) ⊗[(R : Type)]
          (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type)) = 1 := by
  haveI := invertible_gamma_toInvSheaf R L
  exact Module.Invertible.finrank_eq_one _ _

/-! ## ★★★★★★★★★段 B の核 —— 商は捻れである -/

open scoped TensorProduct in
/-- ★★★★★**`Γ(L)` は分数体へ単射である**(＝torsion-free)。

★実測(2026-08-28): `Module.Flat` は `Module.Invertible` から instance で出るが、
`NoZeroSMulDivisors` への橋は instance としては**無い**。
★★代わりに `Module.Flat.tensorProduct_mk_injective` が直接使える。 -/
theorem gamma_toFractionRing_injective (R : CommRingCat.{0}) [IsDomain (R : Type)]
    (L : AInv (Spec R)) :
    Function.Injective
      ((TensorProduct.mk (R : Type) (FractionRing (R : Type))
        (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type)) 1) := by
  haveI := invertible_gamma_toInvSheaf R L
  exact Module.Flat.tensorProduct_mk_injective _ _ _

open scoped TensorProduct in
/-- ★★★★★★★★★**商 `Γ(L)/(R·s)` は捻れである**(`s ≠ 0`)。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★これが台帳 `arakelov-degF-finite-places` の**段 B の核**である。

★機構は 3 つの組み合わせ:

| 道具 | 役割 |
|---|---|
| `finrank_fractionRing_gamma` | 分数体上で rank 1 |
| `gamma_toFractionRing_injective` | `Γ(L) → K ⊗ Γ(L)` が単射 |
| `IsFractionRing.div_surjective` | `c = a/y` に分解 |

★★`s ≠ 0` なら `1 ⊗ s ≠ 0` で、rank 1 なので `1 ⊗ m = c • (1 ⊗ s)`。
`c = a/y` の分母を払って単射性を使えば `y • m = a • s` となる。 -/
theorem exists_smul_mem_span_gamma (R : CommRingCat.{0}) [IsDomain (R : Type)]
    (L : AInv (Spec R))
    (s m : (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type))
    (hs : s ≠ 0) :
    ∃ r ∈ nonZeroDivisors (R : Type), r • m ∈ Submodule.span (R : Type) {s} := by
  have hinj := gamma_toFractionRing_injective R L
  have hts : (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] s ≠ 0 := by
    intro h
    refine hs (hinj ?_)
    show (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] s
      = (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] 0
    rw [h, TensorProduct.tmul_zero]
  obtain ⟨c, hc⟩ := (finrank_eq_one_iff_of_nonzero' _ hts).1 (finrank_fractionRing_gamma R L)
    ((1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] m)
  obtain ⟨a, y, hy, hay⟩ := IsFractionRing.div_surjective (R : Type) c
  refine ⟨y, hy, ?_⟩
  have hy0 : (algebraMap (R : Type) (FractionRing (R : Type))) y ≠ 0 :=
    IsFractionRing.to_map_ne_zero_of_mem_nonZeroDivisors hy
  have h1 : ∀ (r : (R : Type))
      (x : (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type)),
      (algebraMap (R : Type) (FractionRing (R : Type))) r
          • ((1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] x)
        = (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] (r • x) := by
    intro r x
    rw [TensorProduct.smul_tmul', smul_eq_mul, mul_one,
      Algebra.algebraMap_eq_smul_one, TensorProduct.smul_tmul]
  have key : (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] (y • m)
      = (1 : FractionRing (R : Type)) ⊗ₜ[(R : Type)] (a • s) := by
    rw [← h1, ← h1, ← hc, ← hay, smul_smul, mul_comm, div_mul_cancel₀ _ hy0]
  have hym : y • m = a • s := hinj key
  rw [hym]
  exact Submodule.mem_span_singleton.2 ⟨a, rfl⟩

/-! ## ★★★★★★★★段 B の最後の道具 —— 有限生成＋捻れ ⇒ 有限 -/

/-- ★★★★★★★★**有限生成な捻れ Z-加群は有限である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★実測(2026-08-28): mathlib には「有限生成＋捻れ ⇒ 有限」の
**直接の補題は無い**(`exact?` も失敗)。
★代わりに PID 上の構造定理 `Module.equiv_directSum_of_isTorsion` を使う:

    M ≅ ⊕ (Z / (p_i ^ e_i))

★★各因子は `Int.quotientSpanEquivZMod` で `ZMod` になり、
添字は `Fintype` なので `DFinsupp.equivFunOnFintype` で有限積になる。

★★★★これが台帳 `arakelov-degF-finite-places` の**段 B の最後の道具**である。 -/
theorem finite_of_finite_isTorsion_int (M : Type) [AddCommGroup M] [Module.Finite ℤ M]
    (hM : Module.IsTorsion ℤ M) : Finite M := by
  obtain ⟨ι, fι, p, hp, e, ⟨eq⟩⟩ := Module.equiv_directSum_of_isTorsion (R := ℤ) (M := M) hM
  haveI hfin : ∀ i : ι, Finite (ℤ ⧸ (Submodule.span ℤ {p i ^ e i})) := by
    intro i
    have hp0 : (p i ^ e i) ≠ 0 := pow_ne_zero _ (hp i).ne_zero
    haveI : NeZero (p i ^ e i).natAbs := ⟨Int.natAbs_ne_zero.2 hp0⟩
    exact Finite.of_equiv _ (Int.quotientSpanEquivZMod (p i ^ e i)).symm.toEquiv
  haveI : Finite (DirectSum ι (fun i => ℤ ⧸ (Submodule.span ℤ {p i ^ e i}))) :=
    Finite.of_equiv _ (DFinsupp.equivFunOnFintype).symm
  exact Finite.of_equiv _ eq.symm.toEquiv

/-! ## ★★★★★★★共通の消滅元 -/

/-- ★★★★★★★**有限生成な捻れ加群には共通の消滅元がある**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★実測(2026-08-28): mathlib に **∃ r ≠ 0, r • M = 0** の直接の補題は無い
(`exact?` 失敗)。生成元の有限集合の**積**で作る。

★★これが台帳 `arakelov-degF-finite-places` の段 B の**ルート 2**の鍵である
——共通の消滅元 `r` が取れれば `M` は `R/(r)`-加群になり、
`R/(r)` が有限なら `Module.finite_iff_finite` で `M` も有限になる。 -/
theorem exists_common_annihilator (R M : Type) [CommRing R] [AddCommGroup M] [Module R M]
    [Module.Finite R M] (h : ∀ m : M, ∃ r ∈ nonZeroDivisors R, r • m = 0) :
    ∃ r ∈ nonZeroDivisors R, ∀ m : M, r • m = 0 := by
  classical
  obtain ⟨S, hS⟩ := (Module.finite_def.mp inferInstance : (⊤ : Submodule R M).FG)
  choose f hf hf0 using h
  refine ⟨∏ s ∈ S, f s, Submonoid.prod_mem _ (fun s _ => hf s), ?_⟩
  intro m
  have hker : (Submodule.span R (S : Set M)) ≤
      LinearMap.ker ((∏ s ∈ S, f s) • LinearMap.id (R := R) (M := M)) := by
    rw [Submodule.span_le]
    intro x hx
    simp only [SetLike.mem_coe, LinearMap.mem_ker, LinearMap.smul_apply, LinearMap.id_apply]
    rw [← Finset.prod_erase_mul S f (Finset.mem_coe.mp hx), mul_smul, hf0 x, smul_zero]
  rw [hS] at hker
  have hm := hker (Submodule.mem_top : m ∈ (⊤ : Submodule R M))
  simpa using hm

/-! ### ★出典の紐付け(`.src`) -/

def AInv.toInvSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(計量つき算術直線束を層の水準の可逆層と見ること)",
    sectionId := "genell-def-1-1-ii" }

def APicMToPicSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM X → PicSheaf X——台帳の段 A3、写像の水準)",
    sectionId := "genell-def-1-1-ii" }

def APicMToClassGroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM (Spec 𝓞_F) → Cl(F)——deg_F の橋の有限素点側)",
    sectionId := "genell-def-1-1-ii" }

def APicMToPicSheafHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM X →* PicSheaf X が群準同型であること——台帳の段 A3)",
    sectionId := "genell-def-1-1-ii" }

def APicMToClassGroupHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM (Spec 𝓞_F) →* Cl(F)——deg_F の橋の有限素点側、群準同型)",
    sectionId := "genell-def-1-1-ii" }

def sheafifyTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(層化はテンソル積と可換であること)",
    sectionId := "genell-def-1-1-ii" }

def invertible_gamma_toInvSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(局所自明な前層加群から可逆 R-加群 Γ(L,⊤)——台帳の段 A)",
    sectionId := "genell-def-1-1-ii" }

def exists_common_annihilator.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(有限生成な捻れ加群の共通の消滅元——段 B のルート 2)",
    sectionId := "genell-def-1-1-ii" }

def exists_common_annihilator.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.finite_def / Finset.prod_erase_mul / Submonoid.prod_mem"
      (.inMathlib "Module.finite_def") 4,
    .implicitStep
      ("★実測(2026-08-28): mathlib に ∃ r ≠ 0, r • M = 0 の直接の補題は無い" ++
       "(exact? 失敗)。生成元の有限集合の積で作る") 4,
    .implicitStep
      ("★★段 B の残り(2026-08-28 実測): " ++
       "(1) N = Γ(L)/(O_F·s) に本補題を適用して r ≠ 0 を得る、" ++
       "(2) Module.IsTorsionBySet.module で N を O_F/(r)-加群にする、" ++
       "(3) Finite (O_F/(r)) ——**Submodule.finiteQuotientOfFreeOfRankEq** " ++
       "(LinearAlgebra/FreeModule/Finite/Quotient.lean:96)で出るが " ++
       "finrank Z の一致が要る、" ++
       "(4) Module.finite_iff_finite で結論") 4 ]

def finite_of_finite_isTorsion_int.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(有限生成な捻れ Z-加群は有限——段 B の最後の道具)",
    sectionId := "genell-def-1-1-ii" }

def finite_of_finite_isTorsion_int.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.equiv_directSum_of_isTorsion(PID 上の構造定理)"
      (.inMathlib "Module.equiv_directSum_of_isTorsion") 4,
    .citation "[mathlib]" "Int.quotientSpanEquivZMod / DFinsupp.equivFunOnFintype"
      (.inMathlib "Int.quotientSpanEquivZMod") 4,
    .implicitStep
      ("★実測(2026-08-28): mathlib に「有限生成＋捻れ ⇒ 有限」の " ++
       "**直接の補題は無い**(exact? も失敗)。PID 上の構造定理を経由する") 4,
    .implicitStep
      ("★★段 B に残るのは **Z への降下**だけである: " ++
       "(1) Γ(L)/(O_F·s) が Z 上有限生成(O_F が有限生成 Z-加群なので " ++
       "Module.Finite.trans)、(2) Z 上捻れ(exists_smul_mem_span_gamma の r ∈ O_F を " ++
       "そのノルムで Z の元に落とす)") 4 ]

def exists_smul_mem_span_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(商 Γ(L)/(R·s) が捻れであること——台帳の段 B の核)",
    sectionId := "genell-def-1-1-ii" }

def exists_smul_mem_span_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "finrank_fractionRing_gamma(分数体上で rank 1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.finrank_fractionRing_gamma") 4,
    .citation "[mathlib]" "Module.Flat.tensorProduct_mk_injective(平坦なら分数体へ単射)"
      (.inMathlib "Module.Flat.tensorProduct_mk_injective") 4,
    .citation "[mathlib]" "finrank_eq_one_iff_of_nonzero' / IsFractionRing.div_surjective"
      (.inMathlib "finrank_eq_one_iff_of_nonzero'") 4,
    .implicitStep
      ("★実測(2026-08-28): Module.Flat は Module.Invertible から instance で出るが " ++
       "NoZeroSMulDivisors への橋は instance としては無い。" ++
       "代わりに Module.Flat.tensorProduct_mk_injective が直接使える") 4,
    .implicitStep
      ("★★段 B に残るのは 1 つだけ: **有限生成＋捻れの O_F-加群は有限**。" ++
       "O_F は有限生成 Z-加群なので Γ(L)/(O_F·s) もそうであり、" ++
       "捻れなので有限") 4 ]

def finrank_fractionRing_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Γ(L) は分数体上で rank 1——段 B の鍵)",
    sectionId := "genell-def-1-1-ii" }

def finrank_fractionRing_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.Invertible の環を変える底変換(instance)"
      (.inMathlib "Module.Invertible") 4,
    .citation "[mathlib]" "Module.Invertible.finrank_eq_one([Free R M] 前提)"
      (.inMathlib "Module.Invertible.finrank_eq_one") 4,
    .implicitStep
      ("★実測(2026-08-28): finrank_eq_one は [Free R M] を前提とするので " ++
       "Dedekind 域上では直接使えないが、環を変える底変換が " ++
       "instance であり、体上の加群は常に free なので、" ++
       "分数体へ上げれば 1 行で出る") 4,
    .implicitStep
      ("★★残るのは 2 段: (1) Γ(L) が torsion-free" ++
       "(Module.Flat は instance で出るが NoZeroSMulDivisors への橋は " ++
       "mathlib に instance としては無い——実測 2026-08-28)、" ++
       "(2) 有限生成＋捻れの Z-加群は有限") 4 ]

def finite_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Γ(L) が有限生成・射影であること——段 B の前提)",
    sectionId := "genell-def-1-1-ii" }

def finite_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.Invertible から Module.Finite / Module.Projective(インスタンス)"
      (.inMathlib "Module.Invertible") 4,
    .implicitStep
      ("★段 B の残りの道筋(実測 2026-08-28): " ++
       "Γ(L) は有限生成・射影(mathlib のインスタンス)、" ++
       "O_F は有限生成 Z-加群なので Γ(L) も有限生成 Z-加群。" ++
       "★★残るのは「Γ(L)/(O_F·s) が捻れ」であり、" ++
       "それは Γ(L) ⊗ F が 1 次元 F-ベクトル空間で s ≠ 0 がそれを張ることから出る。" ++
       "★★★有限生成＋捻れの Z-加群は有限") 4 ]

def exists_ne_zero_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(算術直線束に非零の大域切断があること——段 B への一歩)",
    sectionId := "genell-def-1-1-ii" }

def exists_ne_zero_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "invertible_gamma_toInvSheaf(Γ(L) が可逆 R-加群＝段 A)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_gamma_toInvSheaf") 4,
    .citation "[mathlib]" "Module.Invertible.linearEquiv(Mv ⊗ M ≃ R)"
      (.inMathlib "Module.Invertible.linearEquiv") 4,
    .implicitStep
      ("★★残っている段 B の核: Γ(L)/(O_F·s) が**有限**であること。" ++
       "★機構は「Γ(L) は有限生成・階数 1 なので商は捻れ、" ++
       "O_F は Dedekind で剰余体が有限だから有限長」である。" ++
       "★★実測(2026-08-28): mathlib に Module.Invertible の API と Ideal.absNorm はあるが、" ++
       "可逆加群から分数イデアルへの橋を自分で繋ぐ必要がある") 4 ]

def APicMToPicSheaf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isLocallyTrivial_sheafify(層化は局所自明性を保つ＝段 A1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_sheafify") 4,
    .citation "[ABC3]" "InvSheaf.ofLocallyTrivial(局所自明な層加群は可逆層)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.InvSheaf.ofLocallyTrivial") 4,
    .citation "[ABC3]" "picSheafEquivClassGroupOF(Pic(Spec 𝓞_F) ≃* Cl(F))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.picSheafEquivClassGroupOF") 4,
    .implicitStep
      ("★測定(2026-08-28): 要る部品が**すべて在庫**だった——" ++
       "isLocallyTrivial_sheafify / InvSheaf.ofLocallyTrivial / PicSheaf.mk_eq_mk / " ++
       "picSheafEquivClassGroupOF。本ファイルはそれらを継ぐだけである") 4,
    .citation "[ABC3]" "sheafifyTensorRight / sheafifyTensorLeft(層化はテンソル積と可換)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.sheafifyTensorRight") 4,
    .implicitStep
      ("★★map_one は要らない——APicM X も PicSheaf X も**群**なので " ++
       "MonoidHom.mk' が map_mul だけから作る") 4,
    .implicitStep
      ("★★★残っている段: アルキメデス成分(計量から無限素点の実数を読む段)と " ++
       "ADiv/APrc ≅ APicM の組み立て(逆向きを含む)が残る") 4 ]

end ABC3.Found.Arakelov
