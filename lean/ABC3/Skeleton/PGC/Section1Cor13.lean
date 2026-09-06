import ABC3.Skeleton.PGC.Section1Defs
import ABC3.Found.PGC.InertiaTransport

/-!
# [pGC] Corollary 1.3 — 惰性群

不分岐判定 `IsUnramifiedAt` と惰性群 `inertia` / `inertiaObject` の**定義**は
`ABC3/Skeleton/PGC/Section1Defs.lean` にある(設計の経緯もそちら)。
本ファイルは主張 `inertia_recoverable` だけを持つ。

## ★2026-09-06(D13 採用): 逸脱の記録

`∀ RD : ResidueCardinality p` をやめ `realResidueCardinality` に固定した。
★理由: `∀ RD` 版は `Check/PGC/Prop12ForallRD.lean` が**原典が偽と述べている命題と同値**だと
証明した(10 例目の退化。「修復が強すぎる主張を作った」型)。
`Interface` に仮説として置いた理由(構成が未構築)は第 1012 で消えている。
★これは**原典の主張を弱めたのではなく、原典が言っていない主張を落とした**ものである。

★本ファイル(Corollary 1.3)は D13 の**下流**である。Cor 1.3 の主語は 2026-09-05 に
既に実物(`residueCardinality` / `subgroupCorrespondence`)へ移してあったので、
statement 自体は今日変えていない。変わったのは**証明が届いたこと**——
`Found/PGC/InertiaTransport.lean::inertia_recoverable_of_prop12` の仮説 `htop`
(「q は Γ_K の同型類だけで決まる」= Proposition 1.2 の q の半分)が、
Proposition 1.2 を実物に固定したことで**無条件に供給できるようになった**
(`Found/PGC/ResidueCardTransport.lean::residueCard_eq_of_absGal_equiv`)。

★これに伴い、定義と主張を 2 ファイルに分けた(理由は `Section1Defs.lean` の冒頭)。
定義そのものは変えていない。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- **[pGC] Corollary 1.3**

原文 (pGC p.3):
> The inertia subgroup I_K ⊆ Γ_K can be determined group-theoretically from Γ_K.

## 条件付き形式化

`RD`(剰余体の元の個数)と `SC`(開部分群と中間体の対応)を仮説に取る。
どちらも原文が §1 で所与として使うが、我々がまだ構成できていないもの。

## 未解決

原文の論拠は Proposition 1.2 を (L, H) に適用することだが、
そこから `I_K` そのものを得る段は書かれていない(上記の設計の経緯を参照)。

## ★★★2026-09-05: **自由な `SC` の下では偽だった**(訂正)——実物に差し替えた

下の 2026-08-14 の監査は「反証には**非正規な** `inertia` が要り、それには
`Γ_K` の非 Galois な3次以上の拡大を構成する必要がある」で止まっていた。
`Found/PGC/QpNonAbelian.lean` でそれを構成した(`ℚ_p(p^{1/ℓ})` は normal でない)
ので、**反証が届いた**:

`H₀ := Gal(K̄/ℚ_p(p^{1/ℓ}))` は開・指数 `ℓ ≥ 2`・**非正規**。`RD` は実物のまま
`SC` だけを `SC.field K H _ := if (HEq H H₀ ∧ H ≠ ⊤) then A else K`
(`A` は剰余体の元の個数が `q^{[Γ:H₀]}` の不分岐拡大)と取ると、判定条件を
満たす開部分群はちょうど `{⊤, H₀}` になり `inertia = H₀`——非正規。だが
`RecoverableFromAbsGal` は内部自己同型でも保たれることを要求するので矛盾。
証明は `Check/PGC/Cor13Degenerate.lean::cor_1_3_statement_false`(`sorry` 無し)。
★`field_top` は満たしたうえで偽になる。

**修理**: `Found/` が実物を構成したので、それについて述べる形に直した——
`residueCardinality p`(`ResidueCardinalityConstruction.lean`)と
`subgroupCorrespondence p`(`SubgroupCorrespondenceConstruction.lean`)。
実物の下では `inertia = Gal(K̄/K^ur)` であり
(`Found/PGC/InertiaIdentification.lean`)、これは**正規**なので反例は塞がる。

★上の設計メモが予告していた「Track B が本物の `ResidueCardinality` を構成した
時点で、ここに依存する全ての statement が一斉に非空虚性の検査を受ける」が、
まさにこの形で起きた。

## 反証可能性(2026-08-14 監査——上のとおり**この結論は覆った**)

**反証できなかった。** それどころか、我々が構成できる唯一の `SC`(`degenerateSC`)の下では
**この定理は証明できる**(`Check.PGC.inertia_recoverable_degenerateSC`、`sorry` 無し)——
`inertia = ⊤` は任意の全射で `⊤` に移るため。Corollary 3.12 では退化 witness が
反例を与えたが、ここでは逆に定理を与える。
反証には**非正規な** `inertia` が要り(`Check.PGC.map_conj_of_normal`)、
それには `Γ_K` の非 Galois な3次以上の拡大を構成する必要がある。
探した範囲は `Check/PGC/RefutationAttempts.lean`。

## ★★★★2026-09-06: 証明が届いた(D13 の下流)

`Found/PGC/InertiaTransport.lean::inertia_recoverable_of_prop12` は
「q が Γ_K の位相群としての同型類だけで決まる」(= Proposition 1.2 の q の半分)を
仮定すれば Corollary 1.3 が従うことを、原文の論拠
"By applying Proposition 1.2 to L and H" のまま形式化していた。
その仮説は `Found/PGC/ResidueCardTransport.lean::residueCard_eq_of_absGal_equiv`
(無条件・`sorry` 無し)がそのまま供給する。
残っていた `.needs` の 2 本(Prop 1.2 の適用・`I_K` の構成)はいずれも
`Found` 側で閉じている(後者は `Section1Defs.lean` の `inertia` そのもの)。
実装は規約(G8)どおり `Found/` にあり
(`Found/PGC/InertiaTransport.lean::inertia_recoverable_real`)、ここはそれへ**委譲する**。 -/
theorem inertia_recoverable :
    (inertiaObject (ABC3.Found.PGC.residueCardinality p)
      (ABC3.Found.PGC.subgroupCorrespondence p)).RecoverableFromAbsGal :=
  ABC3.Found.PGC.inertia_recoverable_real

def inertia_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def inertia_recoverable.needs : List ProofObligation :=
  [ .derivation "Proposition 1.2 を (L, H) に適用して q_L を得る段" 3,
    .implicitStep
      "不分岐判定から惰性群 I_K そのものを構成する段(原文に無い。ここでは開かつ不分岐な部分群の共通部分として構成した)" 3 ]

end ABC3.Skeleton.PGC
