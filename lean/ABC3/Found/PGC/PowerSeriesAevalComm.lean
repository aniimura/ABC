import Mathlib.RingTheory.PowerSeries.PiTopology
import Mathlib.RingTheory.PowerSeries.Substitution

/-!
# 連続な環準同型は `PowerSeries.aeval` と可換(`sorry` 無し)

`PowerSeries.aeval`(冪級数の評価、`trunc`のtruncation-limitとして
定義される)と、係数環`A`を固定する連続な代数準同型`σ:S→ₐ[A]S`が
可換であること——`σ(aeval hz g) = aeval hz' g`(`hz':HasEval (σ z)`)
——を示す、Lubin-Tate理論に一切依存しない純粋な位相環論の一般補題。

## 使い道

`Found/PGC/LubinTateActionInjective.lean`(単射性)・`Found/PGC/
LubinTateActionBijective.lean`(全単射)で確立した内容と組み合わせ、
`Gal(K.closure/K.carrier)`のLubin-Tate作用への同変性
(`σ(a·x)=a·σ(x)`)を確立する際の鍵となる一般補題として使う見込み。

## 証明

各truncation`(trunc N g:Polynomial A)`での可換性は
`Polynomial.aeval_algHom_apply`(多項式の評価と代数準同型の可換性、
有限和なので自明)から直ちに出る。`trunc N g→g`
(`PowerSeries.WithPiTopology.tendsto_trunc_atTop`)+
`PowerSeries.continuous_aeval`+`σ`の連続性を組み合わせ、
`tendsto_nhds_unique`で極限へ持ち上げる——`Found/PGC/
LubinTateIdentityLaw.lean::aeval_formalGroupLaw_eq_of_snd_eq_zero`
と同じtruncation-limit手法の、より単純な(等式のみ、不等式を経由
しない)版。
-/

namespace ABC3.Found.PGC

open scoped PowerSeries.WithPiTopology

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**連続な代数準同型は`PowerSeries.aeval`
と可換**: `σ(aeval hz g) = aeval hz' g`。 -/
theorem algHom_aeval_powerSeries_comm {A S : Type*} [CommRing A] [CommRing S] [Algebra A S]
    [UniformSpace A] [IsUniformAddGroup A] [IsTopologicalSemiring A]
    [UniformSpace S] [IsUniformAddGroup S] [T2Space S] [CompleteSpace S] [IsTopologicalRing S]
    [IsLinearTopology S S] [ContinuousSMul A S]
    (σ : S →ₐ[A] S) (hσcont : Continuous σ)
    (g : PowerSeries A) {z : S} (hz : PowerSeries.HasEval z) (hz' : PowerSeries.HasEval (σ z)) :
    σ (PowerSeries.aeval hz g) = PowerSeries.aeval hz' g := by
  have htrunc : ∀ N : ℕ,
      σ (PowerSeries.aeval hz ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A)) =
      PowerSeries.aeval hz' ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A) := by
    intro N
    rw [PowerSeries.aeval_coe, PowerSeries.aeval_coe, Polynomial.aeval_algHom_apply]
  have h1 : Filter.Tendsto
      (fun N => σ (PowerSeries.aeval hz ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A)))
      Filter.atTop (nhds (σ (PowerSeries.aeval hz g))) :=
    hσcont.continuousAt.tendsto.comp
      ((PowerSeries.continuous_aeval hz).continuousAt.tendsto.comp
        (PowerSeries.WithPiTopology.tendsto_trunc_atTop A g))
  have h2 : Filter.Tendsto
      (fun N => PowerSeries.aeval hz' ((PowerSeries.trunc N g : Polynomial A) : PowerSeries A))
      Filter.atTop (nhds (PowerSeries.aeval hz' g)) :=
    (PowerSeries.continuous_aeval hz').continuousAt.tendsto.comp
      (PowerSeries.WithPiTopology.tendsto_trunc_atTop A g)
  rw [funext htrunc] at h1
  exact tendsto_nhds_unique h1 h2

end ABC3.Found.PGC
