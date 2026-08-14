---
name: abc3-disputed-locus-cor312
description: "公開されている争点は [IUTchIII] Thm 3.11→Cor 3.12 に集中し、その指摘の形が G2/G3 ゲートと一致する。"
metadata: 
  node_type: memory
  type: reference
  originSessionId: fbb27ed5-3b09-4626-a97a-bdf85c79f517
  modified: 2026-08-14T01:53:40.066Z
---

`D:\Math_ABC3\ResearchPaper\0_Source` には二次文献が2本ある——`Scholze-Stix - Why abc is still a conjecture.txt`(693行)と `Joshi - Final Report on the Mochizuki-Scholze-Stix Controversy.txt`(551行)。争点は **[IUTchIII] Theorem 3.11 → Corollary 3.12** に集中している。

指摘の**形**が重要:

- Scholze–Stix(l.631): 「the critical [IUTT-3, Theorem 3.11] does not become false, but **trivial**」——偽ではなく自明になる。
- Joshi(l.28): Arithmetic Holomorphic Structure が「2つ以上あることを定量的に主張するには不十分」——存在の witness が足りない。

**Why:** この2つは、Lean の実演から独立に導いた statement ゲートと同じものである——「偽でなく自明」= **G3**(負の対照)、「witness が足りない」= **G2**(非空虚witness)。形式化が寄与しうるのは元々この種の問いなので、ここが本プロジェクトの最大のレバレッジ点になる。

**How to apply:** Phase 2 の骨格展開は Thm 3.11 / Cor 3.12 への到達を優先する。ただし **二次文献の結論を自分の結論として引用しない**——彼らの指摘は「どこを見るべきか」の道具として使い、判定は自分の機械検査で行う。また、自明化が「原典のせい」か「我々の `Interface` が弱すぎるせい」かは、`Interface` の各フィールドが原典の逐語に裏付けられている(G1)ときにのみ切り分けられる。計画全体は [[abc3-plan-two-track]]。
