import ABC3.Skeleton.PGC.Section1Cor13

/-!
# [pGC] Corollary 1.3 — 退化は消えていない、移動した

`ABC3/Check/PGC/InertiaDegeneracy.lean` で、`I_K` を自由なデータとして置くと
`⊥` でも `⊤` でも Corollary 1.3 が成り立ってしまうことを示した。
`ABC3/Skeleton/PGC/Section1Cor13.lean` では `I_K` を原文の不分岐判定から
**構成**することでこの自由度を消した。

このファイルは、**それでも退化は消えていない**ことを示す——
自由度は `Interface` の `ResidueCardinality` へ**移動した**。

退化した `RD`(`card := fun _ => p`。`isPrimePow` は f = 1 で満たす)を与えると、
不分岐条件は `p = p^[Γ_K:H]` となり、成り立つのは指数 1 のときだけ。
したがって共通部分は `⊤` に潰れる。

## これは設計の失敗ではない

退化の可能性が **1 箇所に局在した**ことが要点。自由なデータのままなら、
退化は各 statement に散らばったまま検査されなかった。
Track B が本物の `ResidueCardinality` を構成した時点で、
それに依存する全ての statement が一斉に非空虚性の検査を受ける。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- 退化した剰余体データ: どの体でも剰余体の元の個数を `p` とする。

`isPrimePow`(q = p^f, f > 0)は f = 1 で満たされるので、
**`Interface` の条件を通過してしまう**。 -/
def degenerateRD : ResidueCardinality p where
  card := fun _ => p
  isPrimePow := fun _ => ⟨1, one_pos, (pow_one p).symm⟩

/-- ★**退化の移動の実証**: 退化した `RD` を与えると惰性群は `⊤` に潰れる。

`I_K` の構成自体は原文どおりなのに、`Interface` の中身が偽物だと
Corollary 1.3 は無内容になる。 -/
theorem degenerateRD_inertia_eq_top (SC : SubgroupCorrespondence p) (K : PAdicLocalField p) :
    inertia (degenerateRD (p := p)) SC K = ⊤ := by
  have hp : Nat.Prime p := Fact.out
  -- 不分岐条件を満たす開部分群は `⊤` しかない
  have hsub : {H : Subgroup K.absGal |
      ∃ hH : IsOpen (H : Set K.absGal), IsUnramifiedAt (degenerateRD (p := p)) SC K H hH}
      ⊆ {⊤} := by
    rintro H ⟨hH, hEq⟩
    -- `IsUnramifiedAt` は `p = p ^ H.index` に潰れる
    have h1 : p = p ^ H.index := by simpa [degenerateRD, IsUnramifiedAt] using hEq
    have h1' : p ^ 1 = p ^ H.index := by rw [pow_one]; exact h1
    have h2 : (1 : ℕ) = H.index := Nat.pow_right_injective hp.two_le h1'
    exact Set.mem_singleton_iff.mpr (Subgroup.index_eq_one.mp h2.symm)
  -- 部分集合の sInf は大きい方以上。`sInf {⊤} = ⊤` なので `⊤ ≤ sInf S`。
  have h4 : (⊤ : Subgroup K.absGal) ≤ inertia (degenerateRD (p := p)) SC K := by
    have h5 := sInf_le_sInf hsub
    rwa [sInf_singleton] at h5
  exact top_le_iff.mp h4

#print axioms degenerateRD_inertia_eq_top

/-!
## 読み

`degenerateRD` は `ResidueCardinality` の条件(`isPrimePow`)を**満たす**。
つまり現在の `Interface` は、本物の剰余体データを一意に決めていない。

決めるのは Track B の仕事である——そこで `Nonempty` の witness が
**構成として**与えられたとき、初めて Corollary 1.3 は内容を持つ。
それまでは「型が付いている」以上のことは言えない。
-/

end ABC3.Check.PGC
