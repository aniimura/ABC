import ABC3.Found.Falt1.AlmostEtale
import Mathlib.RepresentationTheory.Homological.GroupCohomology.LowDegree

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
restriction-corestriction ではなくこの特定の平均化だけ**であり、それは
コサイクルの水準で明示的に書ける(`AlmostEtale.lean` の `transfer_H1`・
`transfer_H2`)。ここではそれを mathlib の `groupCohomology.H1`・`H2`
そのものに関する主張へ翻訳する。

## 残っている部分(正直な記録)

- **一般次数 `i>0`**: `H^1`・`H^2` は本ファイルで閉じた(`Theorem 2.2`・
  `2.3` が使うのはこの2つの次数)。一般の `i` は `inhomogeneousCochains`
  の水準で同じ平均化(`h(g₁..g_{n-1}) := Σ_x (g₁⋯g_{n-1}x)(b)·f(g₁..g_{n-1},x)`)
  を書くか、coinduced 加群と Shapiro の補題(`groupCohomology.coindIso`、
  mathlib にある)経由で一括に処理する必要がある。
- **Galois の仮定から `Σ_{g∈G} g = tr_{B/A}` を出す部分**、および局所化
  `B[1/p]` の水準から `B` へ降ろす部分。
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

end ABC3.Found.Falt1
