/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.QTorsion
import ABC3.Found.GaloisRep.QDomain
import Mathlib.FieldTheory.Galois.Infinite

/-!
# [GenEll] Lemma 3.2, (i) —— **`G_K`-安定な直線の二者択一**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★原文は証明を [FC] に帰しているが、直前の段落に議論が書いてある

`Skeleton/GenEll/Section3.lean` の `lemma_3_2.needs` はこう書いていた:

> ★(i) は `Interface` の `stableLine_dvd_or_cyclotomic` を**そのまま受けている**
> （原文が証明を与えず [FC] Ch. III, Cor. 7.3 に帰しているため）

★★**しかし原文は Lemma 3.2 の直前の段落で議論そのものを書いている**。原文 (GenEll p.15):

> Thus, the above exact sequence splits if and only if qE has an l-th root in K. Note

★★★**Tate 一意化 `E(K̄) ≅ K̄ˣ/q^ℤ`（`Found/GaloisRep/` に構成済み）があれば、
この段は初等的な群論になる。** 本ファイルはそれを書く。

## ★★★★機構——3 行で終わる

`M_l(E) = (K̄ˣ/q^ℤ)[l]` の元 `x` に対し `x^l = q^k`（`k ∈ ℤ`）である
（`QTorsion.lean` の座標 `⟨ζ, π⟩`）。

* `k ≡ 0 (mod l)` の側が **`𝔽_l(1)`**（`μ_l` の方向）である。
* `k ≢ 0 (mod l)` なら、`l` が素数だから Bézout で `k·u + l·w = 1` が取れ、
  `π ≝ x^u·q^w` は **`π^l = q`** を満たす（`hπl`）。しかも `π^k = x` なので、
  `x` の直線が安定なら `π` の直線も安定である（`hstab'`）。

そして**安定性が `π` を `K` へ降ろす**（`sigma_eq_of_stable`）——

    σπ = π^c·q^n  かつ  (σπ)^l = q  ⟹  q^{c+ln-1} = 1  ⟹  c + ln = 1  ⟹  σπ = π

★`q` が 1 の冪根でないことだけを使う。★★`Gal(K̄/K)` のすべてで固定されるので
`π ∈ K`（`InfiniteGalois.mem_range_algebraMap_iff_fixed`——**有限次元を要らない版**）、
したがって `q = π^l` で `l ∣ v_K(q)`。

## ★★★★★これが原文の「either ... or」である

| 原文 | 本ファイル |
|---|---|
| `v_K(q_E) ∈ l·ℤ` | ★`lemma_3_2_i`（`k ≢ 0 mod l` の側） |
| `N = 𝔽_l(1)` | `stable_of_mem_mu`（`k ≡ 0 mod l` の側——★負の対照） |

## ★★逸脱の記録（CLAUDE.md の「逸脱」）

★**`M_l(E)` を `E` の捩れ点としてではなく `K̄ˣ/q^ℤ` の `l`-捩れとして扱う。**
両者を繋ぐのは Tate 一意化（`Found/GaloisRep/TateCurveWitness.lean` の
`uniformization_of_split`）であり、**それは既に建っている**。
★★一意化を `G_K`-同変に持ち上げる段は本ファイルには入っていない
——ここは `K̄ˣ/q^ℤ` の側だけを閉じる。
★★★★**2026-08-27——その段は `Found/GenEll/Lemma32Tate.lean` で閉じた**
（`lemma_3_2_i_tate_all`。一意化は `tatePhiAddEquivAll`、Galois 作用は mathlib の `Point.map`）。

★★★**`L/K` は `IsGalois` として受ける**（原文の `K̄/K`）。
`InfiniteGalois` 版を使うので**有限次元は要らない**。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep

/-! ## ★★★★★★安定性が `π` を `K` へ降ろす -/

/-- ★★★★★**[GenEll] Lemma 3.2, (i) の核** —— `π^l = q` かつ直線が `σ` で安定なら `σπ = π`。

原文 (GenEll p.15):
> Thus, the above exact sequence splits if and only if qE has an l-th root in K. Note

★★機構: `σπ = π^c·q^n` の両辺を `l` 乗すると `q = q^{c+ln}`、
`q` が 1 の冪根でないので `c + ln = 1`、★したがって `σπ = π^{c+ln} = π`。

★★★**`c` が可逆であることすら要らない**——整数 1 本の等式で決まる。 -/
theorem sigma_eq_of_stable {K L : Type} [Field K] [Field L] [Algebra K L]
    {l : ℕ} (q : Kˣ) (π : Lˣ)
    (hπ : π ^ l = Units.map (algebraMap K L : K →* L) q)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    (σ : L ≃ₐ[K] L) (c n : ℤ)
    (hst : Units.map (σ : L →* L) π
      = π ^ c * (Units.map (algebraMap K L : K →* L) q) ^ n) :
    Units.map (σ : L →* L) π = π := by
  set Q : Lˣ := Units.map (algebraMap K L : K →* L) q with hQ
  have hσq : Units.map (σ : L →* L) Q = Q := by
    ext; simp [hQ, AlgEquiv.commutes]
  have h1 : (Units.map (σ : L →* L) π) ^ l = Q := by
    rw [← map_pow, hπ, hσq]
  have hpow : ∀ (a : ℤ), π ^ a * Q ^ n = π ^ (a + l * n) := by
    intro a
    rw [zpow_add, ← hπ, ← zpow_natCast π l, ← zpow_mul]
  have h2 : (π ^ c * Q ^ n) ^ l = Q ^ (c + l * n) := by
    rw [hpow c, ← zpow_natCast _ l, ← zpow_mul, mul_comm _ (l:ℤ), zpow_mul, zpow_natCast, hπ]
  rw [hst, h2] at h1
  have h3 : Q ^ (c + l * n - 1) = 1 := by
    rw [zpow_sub, h1, zpow_one, mul_inv_cancel]
  have hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0 := by
    intro j hj
    refine hqinf j ?_
    have hz : Units.map (algebraMap K L : K →* L) (q ^ j) = 1 := by
      rw [map_zpow]; exact hj
    have hv := congrArg Units.val hz
    simp only [Units.coe_map, Units.val_one, MonoidHom.coe_coe] at hv
    exact Units.ext ((algebraMap K L).injective (by simpa using hv))
  have h4 : c + l * n = 1 := by have := hQinf _ h3; omega
  rw [hst, hpow c, h4, zpow_one]

/-- ★★★★★★**安定な `π` は `K` の元の `l` 乗根である**。

原文 (GenEll p.15):
> Thus, the above exact sequence splits if and only if qE has an l-th root in K. Note

★`sigma_eq_of_stable` を全 `σ` に当てると `π` は `Gal(L/K)` 固定になる。
★★`InfiniteGalois.mem_range_algebraMap_iff_fixed` は**有限次元を要らない**ので
`L = K̄` でそのまま使える。 -/
theorem exists_root_of_stable {K L : Type} [Field K] [Field L] [Algebra K L] [IsGalois K L]
    {l : ℕ} (q : Kˣ) (π : Lˣ)
    (hπ : π ^ l = Units.map (algebraMap K L : K →* L) q)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) π
      = π ^ c * (Units.map (algebraMap K L : K →* L) q) ^ n) :
    ∃ p : Kˣ, p ^ l = q := by
  have hfix : ∀ σ : L ≃ₐ[K] L, σ (π : L) = (π : L) := by
    intro σ
    obtain ⟨c, n, hst⟩ := hstab σ
    exact congrArg Units.val (sigma_eq_of_stable q π hπ hqinf σ c n hst)
  obtain ⟨p0, hp0⟩ := (InfiniteGalois.mem_range_algebraMap_iff_fixed (π : L)).2 hfix
  have hp0ne : p0 ≠ 0 := by
    intro h; apply π.ne_zero; rw [← hp0, h, map_zero]
  refine ⟨Units.mk0 p0 hp0ne, ?_⟩
  have hmap : Units.map (algebraMap K L : K →* L) (Units.mk0 p0 hp0ne) = π :=
    Units.ext (by simpa using hp0)
  have hh : Units.map (algebraMap K L : K →* L) ((Units.mk0 p0 hp0ne) ^ l)
      = Units.map (algebraMap K L : K →* L) q := by rw [map_pow, hmap, hπ]
  exact Units.ext ((algebraMap K L).injective (by
    have := congrArg Units.val hh; simpa using this))

/-! ## ★★★★★★★★`Lemma 3.2, (i)` -/

/-- ★★★★★★★★**[GenEll] Lemma 3.2, (i)** —— `G_K`-安定な直線の二者択一。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★★★**`M_l(E)` の元 `x` は `x^l = q^k` で表される**（`QTorsion.lean` の座標）。
`k ≡ 0 (mod l)` の側が `𝔽_l(1)`（`stable_of_mem_mu`）であり、
★本定理は**もう一方の側**を扱う: `l ∤ k` なら `l ∣ v_K(q_E)`。

★★機構は Bézout `k·u + l·w = 1` で `π ≝ x^u·q^w` に正規化すること。
このとき `π^l = q` かつ `π^k = x` なので直線が保たれ、
`exists_root_of_stable` が効いて `q` が `K` で `l` 乗になる。 -/
theorem lemma_3_2_i {K L : Type} [Field K] [Field L] [Algebra K L] [IsGalois K L]
    {l : ℕ} (hl : l.Prime) (q : Kˣ) (v : Kˣ →* Multiplicative ℤ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    (x : Lˣ) (k : ℤ)
    (hx : x ^ l = (Units.map (algebraMap K L : K →* L) q) ^ k)
    (hk : ¬ ((l:ℤ) ∣ k))
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) x
      = x ^ c * (Units.map (algebraMap K L : K →* L) q) ^ n) :
    (l:ℤ) ∣ vAdd v q := by
  set Q : Lˣ := Units.map (algebraMap K L : K →* L) q with hQ
  have hp : Prime (l : ℤ) := Nat.prime_iff_prime_int.mp hl
  obtain ⟨a, b, hab⟩ := (Prime.coprime_iff_not_dvd hp).mpr hk
  set u : ℤ := b with hu
  set w : ℤ := a with hw
  have hbez : k * u + l * w = 1 := by rw [hu, hw]; linarith [hab]
  set π : Lˣ := x ^ u * Q ^ w with hπdef
  have hxl : ∀ s : ℤ, (x ^ s) ^ l = Q ^ (k * s) := by
    intro s
    rw [← zpow_natCast (x ^ s) l, ← zpow_mul, mul_comm s (l:ℤ), zpow_mul, zpow_natCast, hx,
      ← zpow_mul, mul_comm k s, mul_comm s k]
  have hπl : π ^ l = Q := by
    rw [hπdef, mul_pow, hxl u, ← zpow_natCast (Q ^ w) l, ← zpow_mul, ← zpow_add]
    rw [show k * u + w * (l:ℤ) = 1 by linarith [hbez], zpow_one]
  have hxlw : x ^ (-((l:ℤ) * w)) = Q ^ (-(k * w)) := by
    rw [show -((l:ℤ)*w) = (l:ℤ) * (-w) by ring, zpow_mul, zpow_natCast, hx, ← zpow_mul]
    congr 1; ring
  have hxpi : π ^ k = x := by
    rw [hπdef, mul_zpow, ← zpow_mul, ← zpow_mul]
    rw [show u * k = 1 + (-((l:ℤ)*w)) by linarith [hbez], zpow_add, zpow_one, hxlw]
    rw [mul_assoc, ← zpow_add, show -(k*w) + w*k = 0 by ring, zpow_zero, mul_one]
  have hstab' : ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) π = π ^ c * Q ^ n := by
    intro σ
    obtain ⟨c, n, hst⟩ := hstab σ
    have hσQ : Units.map (σ : L →* L) Q = Q := by ext; simp [hQ, AlgEquiv.commutes]
    have hlhs : Units.map (σ : L →* L) π = x ^ (c * u) * Q ^ (n * u + w) := by
      rw [hπdef, map_mul, map_zpow, map_zpow, hst, hσQ, mul_zpow, ← zpow_mul, ← zpow_mul,
        mul_assoc, ← zpow_add]
    have hkey : π ^ (k * (c * u)) = x ^ (c * u) := by rw [zpow_mul, hxpi]
    exact ⟨k * (c * u), n * u + w, by rw [hkey]; exact hlhs⟩
  obtain ⟨p, hpl⟩ := exists_root_of_stable q π hπl hqinf hstab'
  refine ⟨vAdd v p, ?_⟩
  rw [← hpl]
  show Multiplicative.toAdd (v (p ^ l)) = _
  rw [map_pow]
  simp [vAdd]

/-! ## ★★★★★★G3 負の対照 -/

/-- ★★★★**負の対照** —— **`𝔽_l(1)` の側では結論が出ない**。

`μ_l` の元 `ζ`（`ζ^l = 1`）が生成する直線は、`σζ ∈ ⟨ζ⟩` である限り安定であり、
そのとき `k = 0` なので `l ∣ k` である。★`lemma_3_2_i` の仮定 `¬ l ∣ k` が破れる。

★★★**これが原文の「either ... or」の or の側**であり、
`lemma_3_2_i` が `¬ l ∣ k` を仮定していることが**空虚でない**ことを示す
——実際に両方の場合が起こる。 -/
theorem stable_of_mem_mu {K L : Type} [Field K] [Field L] [Algebra K L]
    {l : ℕ} (q : Kˣ) (ζ : Lˣ) (hζ : ζ ^ l = 1)
    (hfix : ∀ σ : L ≃ₐ[K] L, Units.map (σ : L →* L) ζ ∈ Subgroup.zpowers ζ) :
    ζ ^ l = (Units.map (algebraMap K L : K →* L) q) ^ (0 : ℤ)
      ∧ ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) ζ
          = ζ ^ c * (Units.map (algebraMap K L : K →* L) q) ^ n := by
  refine ⟨by rw [hζ, zpow_zero], fun σ => ?_⟩
  obtain ⟨c, hc⟩ := hfix σ
  exact ⟨c, 0, by rw [zpow_zero, mul_one, ← hc]⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Lemma 3.2` は (i)(ii) の 2 条あり、
また本ファイルは `M_l(E)` を `K̄ˣ/q^ℤ` の側で閉じただけで、
一意化を `G_K`-同変に持ち上げる段は入っていない。
★★★★**その段は 2026-08-27 に `Lemma32Tate.lean` で閉じた**——
本ファイルの条つき `.src` の限定は**本ファイルの範囲の話**であり、項目の残高ではない。 -/

def sigma_eq_of_stable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(安定性が π を K へ降ろす——核)",
    sectionId := "genell-lemma-3-2" }

def exists_root_of_stable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(安定な π は K の元の l 乗根)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(Kˣ/q^ℤ の側——G_K 同変な一意化は含まない)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_root_of_stable(安定な π は K の元の l 乗根)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_root_of_stable") 15,
    .citation "[mathlib]" "InfiniteGalois.mem_range_algebraMap_iff_fixed(固定元は K から来る——有限次元不要)"
      (.inMathlib "InfiniteGalois.mem_range_algebraMap_iff_fixed") 15,
    .citation "[ABC3]" "qQuot_torsion_addEquiv((Kˣ/q^ℤ)[N] ≅ (ℤ/N)² ——座標の出どころ)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.qQuot_torsion_addEquiv") 15,
    .implicitStep
      ("★★★原文は Lemma 3.2 の証明を [FC] Ch. III, Cor. 7.3 に帰しているが、" ++
       "直前の段落(「the above exact sequence splits if and only if qE has an l-th " ++
       "root in K」)に議論そのものが書いてある。Tate 一意化があれば初等的な群論になる") 15,
    .implicitStep
      ("★逸脱: M_l(E) を E の捩れ点としてではなく K̄ˣ/q^ℤ の l-捩れとして扱う。" ++
       "両者を繋ぐ Tate 一意化(uniformization_of_split)は建っているが、" ++
       "それを G_K-同変に持ち上げる段は本ファイルには入っていない。" ++
       "★★★その段は 2026-08-27 に Lemma32Tate.lean で閉じた（lemma_3_2_i_tate_all）") 15,
    .implicitStep
      ("★★機構: Bézout k·u + l·w = 1 で π ≝ x^u·q^w に正規化すると π^l = q かつ π^k = x。" ++
       "そこで σπ = π^c·q^n の両辺を l 乗すると q^{c+ln-1} = 1、q が 1 の冪根でないので" ++
       "c + ln = 1、したがって σπ = π") 15 ]

def stable_of_mem_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(★G3 負の対照——𝔽_l(1) の側では結論が出ない)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
