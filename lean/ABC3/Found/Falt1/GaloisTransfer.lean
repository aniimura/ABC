import ABC3.Found.Falt1.AlmostEtale
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree
import Mathlib.RepresentationTheory.Homological.GroupCohomology.Shapiro
import Mathlib.RepresentationTheory.Homological.GroupCohomology.Functoriality

/-!
# [Falt1] Theorem 2.4(ii) の群コホモロジー側——transfer(平均化)、Found(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §2、物理 p.8
(印字 p.261)。

内容 (Falt1 p.8、260dpi 目視): Suppose B is an almost étale covering of A.
… (ii) Suppose a finite group G operates on B, such that B[1/p] is a Galois
covering of A[1/p] with group G. If M is any B-module with a semilinear
G-action, then m annihilates all higher cohomology H^i(G,M), i>0.

## この節の位置づけ

Faltings 自身の証明は「`m`が`A/tr_{B/A}(B)`を零化することを示せば十分」に
帰着する——これは `Definition 2.1` 直後の remark (iii) そのもので、
`Found/Falt1/AlmostEtale.lean` の `remark_iii_trace_identity`・
`trace_ideal_pow_mem_traceIdeal` で証明済みである。Galois の場合
`tr_{B/A} = Σ_{g∈G} g` なので、remark(iii) は
「**`Σ_{g∈G} g(b) = p^ε` を満たす `b∈B` が存在する**」を与える。

そこから先は古典的な **transfer(平均化)の議論**である。以前のセッションで
「mathlib の `groupCohomology` に一般の transfer 定理(restriction-
corestriction)が無い」ことを確認して壁として報告していたが、**必要なのは
restriction-corestriction ではなくこの特定の平均化だけ**であり、
本ファイルで**全ての次数 `i>0` について証明した**
(`transfer_groupCohomology_smul_eq_zero`)。加えて 2.4(ii) 後半の
「`M^G/tr_G(M)` についても同様」も `transfer_invariants_mem_trace` で
閉じている。

証明の骨格は`AlmostEtale.lean`の`ext_smul_eq_zero_of_almost_split`
(remark 2.1(v)の almost 版で使った「almost split ⟹ almost 消滅」)と
**全く同じ形**である:
1. `Σ_{g∈G} g(b) = c` から、`M` が coinduced 加群
   `Coind_{1}^{G}(M)` の "almost direct summand"(`s≫μ = c•𝟙`)である
   ことを明示的に構成する(`coind_almost_split`——`s(m) := (x ↦ x•m)`、
   `μ(F) := Σ_x x•(b•F(x⁻¹))`、合成が `c•` になるのは semilinear 性
   `g•(b•x) = (g•b)•(g•x)` から)。
2. coinduced 加群のコホモロジーは正次数で消える(Shapiro の補題
   `groupCohomology.coindIso` + 自明群のコホモロジーの消滅
   `isZero_groupCohomology_succ_of_subsingleton`、どちらも mathlib)。
3. 関手性(`groupCohomology.map_id_comp`・`π_map_apply`)で
   `c•` が 0 を経由することを言う。

`H^1`・`H^2` については、より具体的なコサイクル水準の明示公式
(`AlmostEtale.lean` の `transfer_H1`・`transfer_H2`——コバウンダリを
与える元を陽に構成する)も併せて用意してある。

## 残っている部分(正直な記録)

**Galois の仮定から `Σ_{g∈G} g = tr_{B/A}` を出す部分**、および局所化
`B[1/p]` の水準から `B` へ降ろす部分。ここが繋がれば `Theorem 2.4(ii)`
の証明は完成する。
-/

namespace ABC3.Found.Falt1

open groupCohomology

universe u

/-- **`H^1(G,M)` の transfer**(mathlib の `groupCohomology.H1` に関する形)。
`M` が `B`-加群、`G` が `B` にも `M` にも作用して semilinear
(`g•(b•x) = (g•b)•(g•x)`)であり、`c : k` が
`algebraMap k B c = Σ_{g∈G} g(b)` を満たすなら、`c` は `H^1(G,M)` の
全ての元を零化する。

Faltings `Theorem 2.4(ii)` では `k := A`・`c := p^ε`——remark(iii) が
まさにこの形の `b` を与える。証明は `transfer_H1`(コサイクル水準の
明示的な平均化)を `H1π_eq_zero_iff` で翻訳するだけ。 -/
theorem transfer_H1_smul_eq_zero {k B M G : Type u} [CommRing k] [CommRing B] [Algebra k B]
    [Group G] [Fintype G] [MulSemiringAction G B]
    [AddCommGroup M] [Module B M] [Module k M] [IsScalarTower k B M]
    [DistribMulAction G M] [SMulCommClass G k M]
    (hsemi : ∀ (g : G) (b : B) (x : M), g • (b • x) = (g • b) • (g • x))
    (b : B) (c : k) (hc : algebraMap k B c = ∑ g : G, g • b)
    (z : H1 (Rep.ofDistribMulAction k G M)) : c • z = 0 := by
  induction z using groupCohomology.H1_induction_on with
  | h x =>
    have hcoc := isCocycle₁_of_mem_cocycles₁ (k := k) (A := M) (fun g => (x : G → _) g) x.2
    obtain ⟨m, hm⟩ := transfer_H1 hsemi b (fun g => (x : G → _) g)
      (fun g h => (hcoc g h).trans (add_comm _ _))
    have hsm : c • H1π (Rep.ofDistribMulAction k G M) x
        = H1π (Rep.ofDistribMulAction k G M) (c • x) := by rw [map_smul]
    rw [hsm, H1π_eq_zero_iff]
    refine ⟨m, ?_⟩
    ext g
    simp only [d₀₁_hom_apply, Rep.ofDistribMulAction_ρ_apply_apply]
    exact (hm g).symm.trans (by rw [← hc, algebraMap_smul]; rfl)

/-- **`H^2(G,M)` の transfer**(`Theorem 2.2` の障害類が住む次数)。
`transfer_H1_smul_eq_zero` と同じ仮定・同じ結論の `H^2` 版。 -/
theorem transfer_H2_smul_eq_zero {k B M G : Type u} [CommRing k] [CommRing B] [Algebra k B]
    [Group G] [Fintype G] [MulSemiringAction G B]
    [AddCommGroup M] [Module B M] [Module k M] [IsScalarTower k B M]
    [DistribMulAction G M] [SMulCommClass G k M]
    (hsemi : ∀ (g : G) (b : B) (x : M), g • (b • x) = (g • b) • (g • x))
    (b : B) (c : k) (hc : algebraMap k B c = ∑ g : G, g • b)
    (z : H2 (Rep.ofDistribMulAction k G M)) : c • z = 0 := by
  induction z using groupCohomology.H2_induction_on with
  | h x =>
    have hcoc := isCocycle₂_of_mem_cocycles₂ (k := k) (A := M) (fun p => (x : G × G → _) p) x.2
    obtain ⟨y, hy⟩ := transfer_H2 hsemi b (fun p => (x : G × G → _) p) hcoc
    have hsm : c • H2π (Rep.ofDistribMulAction k G M) x
        = H2π (Rep.ofDistribMulAction k G M) (c • x) := by rw [map_smul]
    rw [hsm, H2π_eq_zero_iff]
    refine ⟨y, ?_⟩
    ext p
    simp only [d₁₂_hom_apply, Rep.ofDistribMulAction_ρ_apply_apply]
    exact (hy p.1 p.2).symm.trans (by rw [← hc, algebraMap_smul]; rfl)

/-! ## 一般次数 `i>0`——coinduced 加群への almost split と Shapiro 経由 -/

open CategoryTheory in
/-- **`M` は coinduced 加群 `Coind_1^G(M)` の "almost direct summand"**。
`Σ_{g∈G} g(b) = algebraMap k B c` を満たす `b` があれば、
`s(m) := (x ↦ x•m)`・`μ(F) := Σ_x x•(b•F(x⁻¹))` が `Rep k G` の射になり、
合成が `c•𝟙` になる。`AlmostEtale.lean` の `hochSectionAlmost_comp_lmul`
(remark 2.1(v)の almost 版で `B` が `B⊗_AB` の almost direct summand で
あることを示した部分)と全く同じ形の構成である。 -/
theorem coind_almost_split {k B M G : Type u} [CommRing k] [CommRing B] [Algebra k B]
    [Group G] [Fintype G] [MulSemiringAction G B]
    [AddCommGroup M] [Module B M] [Module k M] [IsScalarTower k B M]
    [DistribMulAction G M] [SMulCommClass G k M]
    (hsemi : ∀ (g : G) (b : B) (x : M), g • (b • x) = (g • b) • (g • x))
    (b : B) (c : k) (hc : algebraMap k B c = ∑ g : G, g • b) :
    ∃ (s : Rep.ofDistribMulAction k G M ⟶ Rep.coind.{u,u,u,u} (⊥ : Subgroup G).subtype
            (Rep.res.{u,u,u,u} (⊥ : Subgroup G).subtype (Rep.ofDistribMulAction k G M)))
      (μ : Rep.coind.{u,u,u,u} (⊥ : Subgroup G).subtype
            (Rep.res.{u,u,u,u} (⊥ : Subgroup G).subtype (Rep.ofDistribMulAction k G M))
          ⟶ Rep.ofDistribMulAction k G M),
      s ≫ μ = c • 𝟙 _ := by
  set ρM := Representation.ofDistribMulAction k G M with hρM
  set φ := (⊥ : Subgroup G).subtype with hφ
  have hbot : ∀ g : (⊥ : Subgroup G), φ g = 1 := fun g => Subgroup.mem_bot.mp g.2
  set ρres := ((Rep.res.{u,u,u,u} φ (Rep.of ρM)).ρ) with hρres
  let sLin : M →ₗ[k] (Representation.coindV φ ρres) :=
    { toFun := fun m => ⟨fun x => x • m, by
        intro g h
        show (φ g * h) • m = ρM (φ g) ((fun x => x • m) h)
        rw [hbot g, one_mul, map_one]
        rfl⟩
      map_add' := by intro x y; ext z; simp [smul_add]
      map_smul' := by intro a x; ext z; simp [smul_comm] }
  let μLin : (Representation.coindV φ ρres) →ₗ[k] M :=
    { toFun := fun F => ∑ x : G, x • (b • (F : G → M) x⁻¹)
      map_add' := by
        intro F F'
        simp only [Submodule.coe_add, Pi.add_apply, smul_add]
        rw [Finset.sum_add_distrib]
      map_smul' := by
        intro a F
        simp only [SetLike.val_smul, Pi.smul_apply, RingHom.id_apply, Finset.smul_sum]
        refine Finset.sum_congr rfl (fun x _ => ?_)
        rw [smul_comm b a ((F : G → M) x⁻¹), smul_comm x a (b • (F : G → M) x⁻¹)] }
  have hsint : ∀ g : G, sLin ∘ₗ ρM g = (Representation.coind φ ρres) g ∘ₗ sLin := by
    intro g
    ext m x
    show x • (g • m) = (x * g) • m
    rw [mul_smul]
  have hμint : ∀ g : G, μLin ∘ₗ (Representation.coind φ ρres) g = ρM g ∘ₗ μLin := by
    intro g
    ext F
    show (∑ x : G, x • (b • (F : G → M) (x⁻¹ * g))) = g • ∑ x : G, x • (b • (F : G → M) x⁻¹)
    rw [Finset.smul_sum]
    refine (Fintype.sum_equiv (Equiv.mulLeft g)
      (fun y : G => (g * y) • (b • (F : G → M) y⁻¹))
      (fun x : G => x • (b • (F : G → M) (x⁻¹ * g))) (fun y => by
        show (g * y) • (b • (F : G → M) y⁻¹) = (g * y) • (b • (F : G → M) ((g * y)⁻¹ * g))
        rw [show ((g * y)⁻¹ * g) = y⁻¹ from by group])).symm.trans ?_
    exact Finset.sum_congr rfl (fun y _ => by rw [mul_smul])
  have hcomp : ∀ m : M, μLin (sLin m) = c • m := by
    intro m
    show (∑ x : G, x • (b • (x⁻¹ • m))) = c • m
    have hterm : ∀ x : G, x • (b • (x⁻¹ • m)) = (x • b) • m := by
      intro x
      rw [hsemi, ← mul_smul, mul_inv_cancel, one_smul]
    rw [Finset.sum_congr rfl (fun x _ => hterm x), ← Finset.sum_smul, ← hc, algebraMap_smul]
  refine ⟨Rep.ofHom ⟨sLin, hsint⟩, Rep.ofHom ⟨μLin, hμint⟩, ?_⟩
  ext m
  exact hcomp m

open CategoryTheory Limits in
/-- **`Theorem 2.4(ii)` 前半の群コホモロジー的内容、全次数版**。
`Σ_{g∈G} g(b) = algebraMap k B c` を満たす `b∈B` があれば、`c` は
**全ての正次数** `H^{n+1}(G,M)` を零化する。Faltings の
「m annihilates all higher cohomology H^i(G,M), i>0」そのもの
(`c := p^ε`、`b` は remark(iii) が与える元)。

証明: `coind_almost_split` で `M` を coinduced 加群の almost direct
summand として捉え、Shapiro の補題(`coindIso`)+ 自明群のコホモロジー
消滅で中間項が消えることを言い、関手性(`map_id_comp`・`π_map_apply`)で
`c•` がその 0 を経由することを示す。`cocyclesMap` が `c•𝟙` を `c•` に
送ることだけは `iCocycles` の単射性(`cyclesMap_i`)で確かめる。 -/
theorem transfer_groupCohomology_smul_eq_zero {k B M G : Type u} [CommRing k] [CommRing B]
    [Algebra k B] [Group G] [Fintype G] [MulSemiringAction G B]
    [AddCommGroup M] [Module B M] [Module k M] [IsScalarTower k B M]
    [DistribMulAction G M] [SMulCommClass G k M]
    (hsemi : ∀ (g : G) (b : B) (x : M), g • (b • x) = (g • b) • (g • x))
    (b : B) (c : k) (hc : algebraMap k B c = ∑ g : G, g • b) (n : ℕ)
    (z : groupCohomology (Rep.ofDistribMulAction k G M) (n+1)) : c • z = 0 := by
  obtain ⟨s, μ, hsμ⟩ := coind_almost_split hsemi b c hc
  have hNzero : IsZero (groupCohomology (Rep.coind.{u,u,u,u} (⊥ : Subgroup G).subtype
      (Rep.res.{u,u,u,u} (⊥ : Subgroup G).subtype (Rep.ofDistribMulAction k G M))) (n+1)) :=
    IsZero.of_iso (isZero_groupCohomology_succ_of_subsingleton _ n)
      (groupCohomology.coindIso _ (n+1))
  have hmapzero : map (A := Rep.ofDistribMulAction k G M) (B := Rep.ofDistribMulAction k G M)
      (MonoidHom.id G) (s ≫ μ) (n+1) = 0 := by
    rw [map_id_comp, hNzero.eq_zero_of_tgt
      (map (A := Rep.ofDistribMulAction k G M) (MonoidHom.id G) s (n+1)), zero_comp]
  induction z using groupCohomology_induction_on with
  | h x =>
    have hcyc : cocyclesMap (A := Rep.ofDistribMulAction k G M) (B := Rep.ofDistribMulAction k G M)
        (MonoidHom.id G) (s ≫ μ) (n+1) x = c • x := by
      rw [hsμ]
      apply (ModuleCat.mono_iff_injective (iCocycles (Rep.ofDistribMulAction k G M) (n+1))).mp
        inferInstance
      have h := HomologicalComplex.cyclesMap_i
        (cochainsMap (MonoidHom.id G)
          (show Rep.res (MonoidHom.id G) (Rep.ofDistribMulAction k G M)
              ⟶ Rep.ofDistribMulAction k G M
            from c • 𝟙 (Rep.ofDistribMulAction k G M))) (n+1)
      have h2 := congrArg (fun (m : _ ⟶ _) => ModuleCat.Hom.hom m x) h
      simp only [ModuleCat.hom_comp, LinearMap.coe_comp, Function.comp_apply] at h2
      rw [h2]
      simp [cochainsMap_id_f_hom_eq_compLeft]
      rfl
    have h1 : map (A := Rep.ofDistribMulAction k G M) (B := Rep.ofDistribMulAction k G M)
        (MonoidHom.id G) (s ≫ μ) (n+1) (π _ (n+1) x)
        = c • π (Rep.ofDistribMulAction k G M) (n+1) x := by
      rw [π_map_apply, hcyc, _root_.map_smul]
    rw [← h1, hmapzero]
    simp

/-- **`Theorem 2.4(ii)` 後半**(「The same holds for `M^G/tr_G(M)`」)。
`G`-不変な `m` について `c•m = tr_G(b•m)` なので、`c` は `M^G/tr_G(M)` を
零化する。semilinear 性から1行で出る。 -/
theorem transfer_invariants_mem_trace {B M G : Type u} [CommRing B]
    [Group G] [Fintype G] [MulSemiringAction G B]
    [AddCommGroup M] [Module B M] [DistribMulAction G M]
    (hsemi : ∀ (g : G) (b : B) (x : M), g • (b • x) = (g • b) • (g • x))
    (b : B) (m : M) (hm : ∀ g : G, g • m = m) :
    (∑ g : G, g • b) • m = ∑ g : G, g • (b • m) := by
  rw [Finset.sum_smul]
  refine Finset.sum_congr rfl (fun g _ => ?_)
  rw [hsemi, hm g]

/-! ## 非空虚性の対照——古典的な「`|G|` が `H^i(G,M)` を零化する」

`B := k`(`G` は自明に作用)・`b := 1` と取ると `Σ_{g∈G} g(b) = |G|` に
なるので、上の2定理はそれぞれ**古典的な「群の位数が正次数のコホモロジーを
零化する」定理**(`|G|` が可逆なら `H^i=0`、という標準的な事実の元になる
主張)に特殊化する。仮定が空虚でないことの確認であると同時に、
主張の向き・符号が正しいことの検算にもなっている。 -/

example {k M G : Type u} [CommRing k] [Group G] [Fintype G] [AddCommGroup M] [Module k M]
    [DistribMulAction G M] [SMulCommClass G k M] (z : H1 (Rep.ofDistribMulAction k G M)) :
    (Fintype.card G : k) • z = 0 := by
  letI : MulSemiringAction G k :=
    { smul := fun _ a => a, one_smul := fun _ => rfl, mul_smul := fun _ _ _ => rfl,
      smul_zero := fun _ => rfl, smul_add := fun _ _ _ => rfl, smul_one := fun _ => rfl,
      smul_mul := fun _ _ _ => rfl }
  haveI := SMulCommClass.symm G k M
  have hsemi : ∀ (g : G) (a : k) (x : M), g • (a • x) = (g • a) • (g • x) := by
    intro g a x
    show g • (a • x) = a • (g • x)
    exact smul_comm g a x
  refine transfer_H1_smul_eq_zero (B := k) hsemi (1 : k) ((Fintype.card G : k)) ?_ z
  show ((Fintype.card G : k)) = ∑ _g : G, (1 : k)
  simp

example {k M G : Type u} [CommRing k] [Group G] [Fintype G] [AddCommGroup M] [Module k M]
    [DistribMulAction G M] [SMulCommClass G k M] (z : H2 (Rep.ofDistribMulAction k G M)) :
    (Fintype.card G : k) • z = 0 := by
  letI : MulSemiringAction G k :=
    { smul := fun _ a => a, one_smul := fun _ => rfl, mul_smul := fun _ _ _ => rfl,
      smul_zero := fun _ => rfl, smul_add := fun _ _ _ => rfl, smul_one := fun _ => rfl,
      smul_mul := fun _ _ _ => rfl }
  haveI := SMulCommClass.symm G k M
  have hsemi : ∀ (g : G) (a : k) (x : M), g • (a • x) = (g • a) • (g • x) := by
    intro g a x
    show g • (a • x) = a • (g • x)
    exact smul_comm g a x
  refine transfer_H2_smul_eq_zero (B := k) hsemi (1 : k) ((Fintype.card G : k)) ?_ z
  show ((Fintype.card G : k)) = ∑ _g : G, (1 : k)
  simp

/-- 全次数版の非空虚性の対照。`|G|` が `H^{n+1}(G,M)` を零化する
(古典的な事実)が `transfer_groupCohomology_smul_eq_zero` から出る。 -/
example {k M G : Type u} [CommRing k] [Group G] [Fintype G] [AddCommGroup M] [Module k M]
    [DistribMulAction G M] [SMulCommClass G k M] (n : ℕ)
    (z : groupCohomology (Rep.ofDistribMulAction k G M) (n+1)) :
    (Fintype.card G : k) • z = 0 := by
  letI : MulSemiringAction G k :=
    { smul := fun _ a => a, one_smul := fun _ => rfl, mul_smul := fun _ _ _ => rfl,
      smul_zero := fun _ => rfl, smul_add := fun _ _ _ => rfl, smul_one := fun _ => rfl,
      smul_mul := fun _ _ _ => rfl }
  haveI := SMulCommClass.symm G k M
  have hsemi : ∀ (g : G) (a : k) (x : M), g • (a • x) = (g • a) • (g • x) := by
    intro g a x
    show g • (a • x) = a • (g • x)
    exact smul_comm g a x
  refine transfer_groupCohomology_smul_eq_zero (B := k) hsemi (1 : k)
    ((Fintype.card G : k)) ?_ n z
  show ((Fintype.card G : k)) = ∑ _g : G, (1 : k)
  simp

end ABC3.Found.Falt1
