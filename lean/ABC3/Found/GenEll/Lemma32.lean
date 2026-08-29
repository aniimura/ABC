/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma32Uniformized
import ABC3.Found.GenEll.Lemma32QuotMu
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Lemma 3.2`（`Found`、項目まるごと）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Lemma 3.2`（Local Rank One Subgroups of `l`-Torsion）を、
**Tate 一意化の言葉で 1 本にまとめる**:

* **(i)** `M_l(E)` の `G_K`-安定な直線は、`𝔽_l(1)`（＝ `μ_l` の像）であるか、
  さもなくば **`l ∣ v_K(q_E)`**。
* **(ii)** `E/μ_l` は**母数 `q_E^l` の Tate 曲線**であり、
  `v_K(q_E^l) = l·v_K(q_E)`——すなわち **`deg_∞(E′) = l·deg_∞(E)`**。
* ★**空虚封じ**: `μ_l` は `E` に**埋め込まれる**（`q` が 1 の冪根でないことが効く）。

## ★★2026-08-29 の前進（明示）—— なぜ今なら項目まるごと取れるのか

`Lemma32QuotMu.lean`（2026-08-27 時点）は `.src` にこう書いていた:

> ★★**項目全体の `.src` は置かない。** …
> **一意化を `G_K`-同変に、また商と両立するように持ち上げる段**が残っている。

★★**前者は 2026-08-27 に解消した**——`Found/GaloisRep/` の 8 段
（`tatePhiAddEquivAll`・`tatePhi_map`）が `G_K`-同変な一意化を**構成**し、
`Lemma32Uniformized.lean` がそれを (i) に接続した。

★★★**後者は逸脱として記録する**（下記）。消費側（`Lemma 3.5`）が
`(ii)` から使うのは **`deg_∞(E′) = l·deg_∞(E)` だけ**であり、
`deg_∞` は `q` の付値で定義されるので、読み替えは消費側に影響しない。

## ★★★★★★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| `E` | `Spec(𝒪_K)` 上の半アーベルスキーム（特殊ファイバーが `(G_m)_k`） | ★**同変な一意化 `Φ : L̄ˣ/q^ℤ ≃+ E` を持つ加法群** | ★原文も `Lemma 3.2` の**前に**この設定を置き、`M_l(E)` の完全列を [FC] Ch. III, Cor. 7.3 に「well-known」として帰している。★★同変一意化は本プロジェクトが**構成済み**（2026-08-27）なので、受けているのは「`E` がその形をしている」ことだけである |
| `μ_l ⊆ E` は**有限平坦部分群スキーム**、`E′ = E/μ_l` は**スキームの商** | スキーム | ★**一意化の側の群の商**（`(L̄ˣ/q^ℤ)/μ_l`） | ★有限平坦群スキームによる商は mathlib にも ABC3 にも無い。★★消費側（`Lemma 3.5`）が使うのは `deg_∞(E′) = l·deg_∞(E)` だけであり、`deg_∞ = v_K(q)·log #k` は付値で決まるので影響しない |
| `N ⊆ M_l(E)` が `𝔽_l(1)` であること | 完全列 `0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0` の部分加群 | ★**`P` が `μ_l` の像に入ること**（`∃ y, y^l = 1 ∧ Φ(mk y) = P`） | ★一意化の下で `𝔽_l(1)` はちょうど `μ_l` の像である。★★`dvd_iff_mem_mu` がこれを `l ∣ k` と**同値**にしており、数値条件ではなく**原文の語**で述べられている |
| `K̄` が代数閉であること | 設定 | ★`hsurj`（`x ↦ x^l` が `L̄ˣ` 上全射） | ★代数閉体では `l` 乗写像は全射である。使うのはこの帰結だけ |

★★★★★★**空虚ではない**——(i) は `Found/GaloisRep/` の同変一意化を実際に使い、
(ii) は核の等式・第 3 同型定理・準同型定理から**構成**されている。
★`mu_inj` が「`μ_l` の像が潰れていない」ことを保証する（潰れていれば
`E/μ_l = E` で `q_{E′} = q_E` となり (ii) が空虚になる）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep Subgroup QuotientGroup

/-! ## ★★★★★`𝔽_l(1)` であることの言い換え -/

/-- ★★★★★**`l ∣ k` ⟺ `P` は `𝔽_l(1)`（＝ `μ_l` の像）に入る**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★`x^l = Q^k` の `k` は代表 `x` の取り方で `l·ℤ` だけ動くので、`k mod l` は
`P` だけで決まる。★★`l ∣ k` なら `y ≝ x·Q^{-k/l}` が `y^l = 1` を満たす。
★★★逆は `q` が 1 の冪根でないことによる。

★★★★これで `Lemma 3.2, (i)` の二者択一を**原文の語（`N = 𝔽_l(1)`）**で書ける。 -/
theorem dvd_iff_mem_mu {L : Type} [Field L] {l : ℕ} (Q : Lˣ)
    (hQinf : ∀ j : ℤ, Q ^ j = 1 → j = 0)
    {E : Type} [AddCommGroup E]
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers Q) ≃+ E)
    (P : E) (x : Lˣ) (k : ℤ)
    (hxP : Φ (Additive.ofMul (QuotientGroup.mk x)) = P)
    (hxk : x ^ l = Q ^ k) :
    (l : ℤ) ∣ k ↔ ∃ y : Lˣ, y ^ l = 1 ∧ Φ (Additive.ofMul (QuotientGroup.mk y)) = P := by
  constructor
  · rintro ⟨m, rfl⟩
    refine ⟨x * (Q ^ m)⁻¹, ?_, ?_⟩
    · rw [mul_pow, hxk, ← zpow_natCast ((Q ^ m)⁻¹) l, inv_zpow, ← zpow_mul, ← zpow_neg,
        ← zpow_add]
      simp [mul_comm]
    · rw [← hxP]
      congr 2
      refine QuotientGroup.eq.2 ?_
      have h : (x * (Q ^ m)⁻¹)⁻¹ * x = Q ^ m := by group
      rw [h]
      exact Subgroup.zpow_mem_zpowers Q m
  · rintro ⟨y, hyl, hyP⟩
    have hmk : (QuotientGroup.mk y : Lˣ ⧸ Subgroup.zpowers Q) = QuotientGroup.mk x :=
      Additive.ofMul.injective (Φ.injective (hyP.trans hxP.symm))
    have hmem : y⁻¹ * x ∈ Subgroup.zpowers Q := QuotientGroup.eq.1 hmk
    obtain ⟨m, hm⟩ := hmem
    have hm' : Q ^ m = y⁻¹ * x := hm
    have hx : x = y * Q ^ m := by rw [hm']; group
    have hQk : Q ^ k = Q ^ ((l : ℤ) * m) := by
      rw [← hxk, hx, mul_pow, hyl, one_mul, ← zpow_natCast (Q ^ m) l, ← zpow_mul, mul_comm]
    have hz : Q ^ (k - (l : ℤ) * m) = 1 := by
      rw [zpow_sub, hQk, mul_inv_cancel]
    have h0 := hQinf _ hz
    exact ⟨m, by omega⟩

/-- ★`q` が `K` で 1 の冪根でなければ、`L` へ送っても 1 の冪根でない。 -/
theorem map_zpow_eq_one_iff {K L : Type} [Field K] [Field L] [Algebra K L] (q : Kˣ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) :
    ∀ j : ℤ, (Units.map (algebraMap K L : K →* L) q) ^ j = 1 → j = 0 := by
  intro j hj
  refine hqinf j ?_
  have h1 : (Units.map (algebraMap K L : K →* L)) (q ^ j) = 1 := by
    rw [map_zpow]; exact hj
  have h2 : (algebraMap K L) ((q ^ j : Kˣ) : K) = 1 := by
    have := congrArg (fun u : Lˣ => (u : L)) h1
    simpa using this
  refine Units.ext ?_
  apply (algebraMap K L).injective
  simpa using h2

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Lemma 3.2` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Lemma 3.2**
(Local Rank One Subgroups of `l`-Torsion)。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

**(i)** `M_l(E)` の `G_K`-安定な直線 `⟨P⟩` は、`𝔽_l(1)`（`μ_l` の像）であるか、
さもなくば **`l ∣ v_K(q_E)`**。

**(ii)** `E/μ_l` は**母数 `q_E^l` の Tate 曲線**であり、
`v_K(q_E^l) = l·v_K(q_E)`——すなわち **`deg_∞(E′) = l·deg_∞(E)`**。

★**空虚封じ**: `μ_l` は `E` に埋め込まれる。

★★仮定の出どころ:

| 仮定 | 原文 |
|---|---|
| `Φ`・`hequiv` | 「`E` は特殊ファイバーが `(G_m)_k` の半アーベルスキーム」——★同変一意化は `Found/GaloisRep/` が**構成済み**（2026-08-27） |
| `hqinf` | `q_E ∈ m_K`（1 の冪根ではない） |
| `hsurj` | `K̄` が代数閉 |
| `hl` | 「`l` be a prime number」 |

★★★逸脱（`E/μ_l` を一意化の側の群の商で読むこと 等）はファイル冒頭の表に記録した。 -/
theorem lemma_3_2 {K L : Type} [Field K] [Field L] [Algebra K L] [IsGalois K L]
    {l : ℕ} (hl : l.Prime) (q : Kˣ) (v : Kˣ →* Multiplicative ℤ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    (hsurj : Function.Surjective (powMonoidHom l : Lˣ →* Lˣ))
    {E : Type} [AddCommGroup E]
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q)) ≃+ E)
    (act : (L ≃ₐ[K] L) → E →+ E)
    (hequiv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ),
      act σ (Φ (Additive.ofMul (QuotientGroup.mk u)))
        = Φ (Additive.ofMul (QuotientGroup.mk (Units.map (σ : L →* L) u)))) :
    (∀ (P : E) (x : Lˣ) (k : ℤ),
        Φ (Additive.ofMul (QuotientGroup.mk x)) = P →
        x ^ l = (Units.map (algebraMap K L : K →* L) q) ^ k →
        (∀ σ : L ≃ₐ[K] L, ∃ c : ℤ, act σ P = c • P) →
        (l : ℤ) ∣ vAdd v q
          ∨ ∃ y : Lˣ, y ^ l = 1 ∧ Φ (Additive.ofMul (QuotientGroup.mk y)) = P)
  ∧ Nonempty (((Lˣ ⧸ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q)) ⧸
        ((powMonoidHom l : Lˣ →* Lˣ).ker.map
          (QuotientGroup.mk' (Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q)))))
      ≃* (Lˣ ⧸ Subgroup.zpowers ((Units.map (algebraMap K L : K →* L) q) ^ l)))
  ∧ vAdd v (q ^ l) = l * vAdd v q
  ∧ (∀ x ∈ (powMonoidHom l : Lˣ →* Lˣ).ker,
        (QuotientGroup.mk' (Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q))) x = 1
          → x = 1) := by
  have hQinf := map_zpow_eq_one_iff (K := K) (L := L) q hqinf
  refine ⟨?_, ⟨quotMuEquiv l _ hsurj⟩, vAdd_pow v l q, mu_inj hl.pos _ hQinf⟩
  intro P x k hxP hxk hstab
  by_cases hk : (l : ℤ) ∣ k
  · exact Or.inr ((dvd_iff_mem_mu _ hQinf Φ P x k hxP hxk).1 hk)
  · exact Or.inl (lemma_3_2_i_of_uniformization hl q v hqinf Φ act hequiv P hstab x k hxP hxk hk)

/-! ## ★出典の紐付け(`.src`)——★★**項目まるごと** -/

def dvd_iff_mem_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(𝔽_l(1) であることの言い換え——l ∣ k ⟺ μ_l の像)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15, item := "Lemma 3.2",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_2_i_of_uniformization((i)——同変一意化を入力として)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_2_i_of_uniformization") 4,
    .citation "[ABC3]" "quotMuEquiv((ii)——E/μ_l は母数 q^l の Tate 曲線)"
      (.inProject "ABC3" "ABC3.Found.GenEll.quotMuEquiv") 4,
    .citation "[ABC3]" "vAdd_pow(v_K(q^l) = l·v_K(q)——deg_∞(E′) = l·deg_∞(E) の実質)"
      (.inProject "ABC3" "ABC3.Found.GenEll.vAdd_pow") 2,
    .citation "[ABC3]" "mu_inj(μ_l は E に埋め込まれる——空虚封じ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mu_inj") 3,
    .otherPaper "[FC]"
      ("Chapter III, Corollary 7.3(M_l(E) の完全列 0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0)。" ++
       "★原文はこれを『well-known』として引く。★★本形式化では一意化 Φ の下で " ++
       "𝔽_l(1) を『μ_l の像』として表しており(dvd_iff_mem_mu)、完全列そのものは要らない") 5,
    .implicitStep
      ("★★★★★★★★逸脱(2026-08-29): 原文は μ_l ⊆ E を**有限平坦部分群スキーム**として、" ++
       "E′ = E/μ_l を**スキームの商**として作る。本形式化は" ++
       "**一意化の側の群の商**((L̄ˣ/q^ℤ)/μ_l)で語る。" ++
       "★有限平坦群スキームによる商は mathlib にも ABC3 にも無い。" ++
       "★★★消費側(Lemma 3.5)が (ii) から使うのは deg_∞(E′) = l·deg_∞(E) だけであり、" ++
       "deg_∞ = v_K(q)·log #k は付値で決まるので、この読み替えは消費側に影響しない") 8,
    .implicitStep
      ("★★★★★2026-08-27 に解消した段: Lemma32QuotMu.lean は" ++
       "『一意化を G_K-同変に…持ち上げる段が残っている』として項目全体の .src を" ++
       "置いていなかった。★その同変一意化は Found/GaloisRep/ の 8 段" ++
       "(tatePhiAddEquivAll・tatePhi_map)で**構成された**。" ++
       "★★本ファイルはそれを (i) に接続し、(ii) と合わせて項目をまるごと取る") 8,
    .implicitStep
      ("★★空虚封じ: mu_inj が『μ_l の像が潰れていない』ことを保証する。" ++
       "潰れていれば E/μ_l = E で q_{E′} = q_E となり (ii) が空虚になる" ++
       "——ここで q が 1 の冪根でないことが効いている") 6 ]

end ABC3.Found.GenEll
